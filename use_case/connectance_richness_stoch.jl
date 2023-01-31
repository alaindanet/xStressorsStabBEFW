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
addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

import Pkg
using DifferentialEquations, BEFWM2, Distributions, ProgressMeter, SparseArrays, LinearAlgebra, DataFrames, CSV
using StatsBase
include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/interaction_strength.jl")

import Random.seed!

seed!(22)

ti = simCS(.4, 90;
           Z = 100,
           d = 0.1, σₑ = .5, ρ = 0.0,
           h = 2.0, c = 1.0,
           r = 1.0, K = 1.0, alpha_ij = 0.5,
           max = 50, last = 10, dt = 0.1,
           fun = stoch_d_dBdt!,
           K_alpha_corrected = true,
           return_sol = false
          )

# Max
FoodWeb(nichemodel, 5, C = .32, tol = .01)
# Min
FoodWeb(nichemodel, 5, C = .16, tol = .01)

# Max
FoodWeb(nichemodel, 10, C = .31, tol = .01)
# Min
FoodWeb(nichemodel, 10, C = .09, tol = .01)

# Max
FoodWeb(nichemodel, 20, C = .24, tol = .01)
# Min
FoodWeb(nichemodel, 20, C = .05, tol = .01)

# Max
FoodWeb(nichemodel, 40, C = .15, tol = .01)
# Min
FoodWeb(nichemodel, 40, C = .03, tol = .01)

# Max
FoodWeb(nichemodel, 80, C = .09, tol = .01)
# Min
FoodWeb(nichemodel, 80, C = .02, tol = .01)

# Max
FoodWeb(nichemodel, 100, C = .07, tol = .01)
# Min
FoodWeb(nichemodel, 100, C = .02, tol = .01)

plot([5, 10, 20, 40, 80, 100], [.32, .31, .24, .15, .09, .07],)

function ctS(S = 10)
    out = (
           min = round((S - 1)/ S^2, digits = 2),
           max = round(0.30 + (.07 - .31) / (100 - 10) * S, digits = 2)
   )
end
FoodWeb(nichemodel, 80, C = .02, tol = .01)
#

# Parameter product
#
#
#
rep = 1:30
S = [5, 10, 20, 40, 60]
C = 0.02:.05:.32
sigma = 0.5
Z = [10, 100, 1000]
ρ = [0, .5, 1]
names = (:rep, :richness, :connectance, :sigma, :rho)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(rep, S, C, sigma, ρ))[:]


# Filter impossible combination of C/S
limitCS = (
           S = S,
           Cmin = round.([(i - 1)/ i^2 + .01 for i in S], digits = 2),
           Cmax = [.32, .31, .24, .15, .09, .07]
          )
goodCSparam_v = [
                 findall(x ->
                         (x.connectance >= limitCS.Cmin[i] && x.connectance <= limitCS.Cmax[i]) && x.richness == limitCS.S[i], param)
                 for i in 1:length(limitCS.S)
                ]
goodCSparam_idxs = reduce(vcat, goodCSparam_v)

bad_param = param[1:length(param) .∉ Ref(goodCSparam_idxs)]
# Select good parameter combination
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


### TEST
#
sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, ),
                               simCS(p.connectance, p.richness;
                                     Z = 100,
                                     d = 0.1, σₑ = p.sigma, ρ = p.rho,
                                     h = 2.0, c = 1.0,
                                     r = 1.0, K = 1.0, alpha_ij = 0.5,
                                     max = 5000, last = 100, dt = 0.1,
                                     fun = stoch_d_dBdt!,
                                     K_alpha_corrected = true,
                                     return_sol = false,
                                     gc_thre = 1
                                    )
                              ),
                         good_param[10 - 5:10]
                        )

sim = @showprogress pmap(p -> merge(
                                    (rep = p.rep, ),
                                    simCS(p.connectance, p.richness;
                                          Z = 100,
                                          d = 0.1, σₑ = p.sigma, ρ = p.rho,
                                          h = 2.0, c = 1.0,
                                          r = 1.0, K = 1.0, alpha_ij = 0.5,
                                          max = 5000, last = 100, dt = 0.1,
                                          fun = stoch_d_dBdt!,
                                          K_alpha_corrected = true,
                                          return_sol = false,
                                          gc_thre = .02
                                         )
                                   ),
                         good_param,
                         batch_size = 100
                        )

df = DataFrame(sim)

serialize("cs_sim.dat", df)
CSV.write("cs_sim.csv", df)
