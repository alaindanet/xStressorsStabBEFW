println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
#flag = "--project=."
println("Workers run with flag: $(flag)")
#addprocs(ncpu - 2, exeflags=flag)
#addprocs(5, exeflags=flag)
#println("Using $(ncpu -2) cores")

import Pkg
using DifferentialEquations, EcologicalNetworksDynamics, Distributions, ProgressMeter
using SparseArrays, LinearAlgebra, DataFrames, CSV
using StatsBase, Arrow
include("../src/stochastic_mortality_model.jl")
include("../src/sim.jl")
include("../src/interaction_strength.jl")

param = DataFrame(Arrow.Table("scripts/param_comb_connectance_richness.arrow"))
param = NamedTuple.(eachrow(param))

println("From $(ARGS[1]) to $(ARGS[2])")

first_sim = parse(Int, ARGS[1])
last_sim = parse(Int, ARGS[2])

if last_sim > size(param, 1)
    last_sim = size(param, 1)
end

println("Running param sim from lines $first_sim to $last_sim")

pm = sample(param)
simCS(pm.connectance, pm.richness;
      Z = 100,
      d = 0.1, σₑ = pm.sigma, ρ = pm.rho,
      h = 2.0, c = 1.0,
      r = pm.enrich.r, K = pm.enrich.K, alpha_ij = 0.5,
      max = 50, last = 10, dt = 0.1,
      fun = stoch_d_dBdt!,
      K_alpha_corrected = true,
      return_sol = false
     )


### TEST
#
sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, productivity = p.enrich.name),
                               simCS(p.connectance, p.richness;
                                     Z = p.Z,
                                     d = 0.1, σₑ = p.sigma, ρ = p.rho,
                                     h = 2.0, c = 1.0,
                                     r = pm.enrich.r, K = pm.enrich.K,
                                     alpha_ij = 0.5,
                                     max = 5000, last = 100, dt = 0.1,
                                     fun = stoch_d_dBdt!,
                                     K_alpha_corrected = true,
                                     return_sol = false,
                                     gc_thre = .2
                                    )
                              ),
                         param[first_sim:last_sim]
                        )

df = DataFrame(sim)

file = string("simCS_test_", first_sim, "_", last_sim, ".arrow")
Arrow.write(file, df)
