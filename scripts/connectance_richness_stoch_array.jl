println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

first_sim = parse(Int, ARGS[1])
last_sim = parse(Int, ARGS[2])

println("Running parameters combination from $(ARGS[1]) to $(ARGS[2])")

ncpu = maximum([length(Sys.cpu_info()), 15])

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -1) cores")


@everywhere import Pkg, Random.seed!
@everywhere using DifferentialEquations, EcologicalNetworksDynamics, SparseArrays
@everywhere using LinearAlgebra, DataFrames
@everywhere using Distributions, ProgressMeter
@everywhere using StatsBase, CSV, Arrow
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../src/interaction_strength.jl")


param = DataFrame(Arrow.Table("/home/bi1ahd/xStressorsStabBEFW/scripts/param_comb_connectance_richness.arrow"))
param = NamedTuple.(eachrow(param))


pm = sample(param)
println("Running warmup: r = $(pm.enrich.r), K = $(pm.enrich.K), σₑ = $(pm.sigma), ρ = $(pm.rho)")

warmup = simCS(pm.connectance, pm.richness;
      Z = 100,
      d = 0.1, σₑ = pm.sigma, ρ = pm.rho,
      h = 2.0, c = 0.0,
      r = pm.enrich.r, K = pm.enrich.K, alpha_ij = 0.5,
      max = 50, last = 10, dt = 0.1,
      K_alpha_corrected = true,
      return_sol = false
     )
println("$(warmup)")


if last_sim > size(param, 1)
    last_sim = size(param, 1)
end
println("Running param sim from lines $first_sim to $last_sim")

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, productivity = p.enrich.name),
                               simCS(p.connectance, p.richness;
                                     Z = p.Z,
                                     d = 0.1, σₑ = p.sigma, ρ = p.rho,
                                     h = 2.0, c = 0.0,
                                     r = p.enrich.r, K = p.enrich.K,
                                     alpha_ij = 0.5,
                                     max = 5000, last = 100, dt = 0.1,
                                     K_alpha_corrected = true,
                                     return_sol = false,
                                     gc_thre = .2
                                    )
                              ),
                         param[first_sim:last_sim],
                         batch_size = 100
                        )

df = DataFrame(sim)
println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

file = string("simCSZ_", first_sim, "_", last_sim, ".arrow")
Arrow.write(file, df)
