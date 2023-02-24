println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

@everywhere import Pkg
@everywhere using DifferentialEquations, BEFWM2, SparseArrays,LinearAlgebra, DataFrames
@everywhere using qistributions, ProgressMeter
@everywhere using StatsBase, CSV, Arrow
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../src/interaction_strength.jl")

import Random.seed!

seed!(22)

rep = 1:20
S = [5, 10, 20, 40, 60]
C = 0.02:.05:.32
sigma = [0.5, 1.0]
Z = [50, 200, 500]
ρ = [0, .5, 1]
names = (:rep, :richness, :connectance, :Z, :sigma, :rho)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(rep, S, C, Z, sigma, ρ))[:]


# Filter impossible combination of C/S
limitCS = (
           S = S,
           Cmin = round.([(i - 1)/ i^2 + .01 for i in S], digits = 2),
           Cmax = [.32, .31, .24, .15, .09, .07]
          )
# Select good combinations
goodCSparam_v = [
                 findall(x ->
                         (x.connectance >= limitCS.Cmin[i] && x.connectance <= limitCS.Cmax[i]) && x.richness == limitCS.S[i], param)
                 for i in 1:length(limitCS.S)
                ]
goodCSparam_idxs = reduce(vcat, goodCSparam_v)
bad_param = param[1:length(param) .∉ Ref(goodCSparam_idxs)]
good_param = param[StatsBase.sample(goodCSparam_idxs, length(goodCSparam_idxs), replace=false)]

# Warm-up
test_i = 2400
simCS(bad_param[test_i].connectance, bad_param[test_i].richness;
      Z = 100,
      d = 0.1, σₑ = bad_param[test_i].sigma, ρ = bad_param[test_i].rho,
      h = 2.0, c = 1.0,
      r = 1.0, K = 1.0, alpha_ij = 0.5,
      max = 50, last = 10, dt = 0.1,
      fun = stoch_d_dBdt!,
      K_alpha_corrected = true,
      return_sol = false
     )

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, ),
                               simCS(p.connectance, p.richness;
                                     Z = p.Z,
                                     d = 0.1, σₑ = p.sigma, ρ = p.rho,
                                     h = 2.0, c = 1.0,
                                     r = 1.8, K = 20.0,
                                     alpha_ij = 0.5,
                                     max = 5000, last = 100, dt = 0.1,
                                     fun = stoch_d_dBdt!,
                                     K_alpha_corrected = true,
                                     return_sol = false,
                                     gc_thre = .2
                                    )
                              ),
                         good_param,
                         batch_size = 100
                        )

df = DataFrame(sim)
println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

Arrow.write("sim_csZ2.arrow", df)
