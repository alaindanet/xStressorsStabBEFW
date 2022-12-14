println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())
println("Using $(ncpu -2) cores")

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 2, exeflags=flag)

@everywhere import Pkg
@everywhere using DifferentialEquations, BEFWM2, Distributions, ProgressMeter, SparseArrays, LinearAlgebra, DataFrames, CSV
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../interaction_strength.jl")
# @everywhere include("src/stochastic_mortality_model.jl")
# @everywhere include("src/sim.jl")
# @everywhere include("src/interaction_strength.jl")

import Random.seed!

seed!(22)

ti = simCS(.1, 10, 100, 2.0, 1.0, 0.5, 5; max = 50, last = 1000, dt = 0.1, return_sol = false)

# Parameter product
#
#
rep = 1:50
S = 10:10:50
C = 0.1:.1:.40
sigma = 0.5:1
names = (:rep, :richness, :connectance, :sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(rep, S, C, sigma))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, ),
                     #simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
                     simCS(p.connectance, p.richness, 100, 2.0, 1.0, p.sigma, 5; max = 5000, last = 1000, dt = 0.1, return_sol = false)
                    ),
           param[1:4]
          )

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, ),
                     #simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
                     simCS(p.connectance, p.richness, 100, 2.0, 1.0, p.sigma, 5; max = 5000, last = 1000, dt = 0.1, return_sol = false)
                    ),
           param,
           batch_size = 100
          )

df = DataFrame(sim)

serialize("cs_sim.dat", df)
CSV.write("cs_sim.csv", df)
