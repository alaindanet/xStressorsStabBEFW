println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using BEFWM2, LinearAlgebra, DifferentialEquations, DataFrames, CSV, Distributed, Distributions, ProgressMeter, BenchmarkTools

ncpu = length(Sys.cpu_info())
println("Using $(ncpu -2) cores")

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
println("Workers run with flag: $(flag)")
addprocs(ncpu - 2, exeflags=flag)

@everywhere import Pkg
@everywhere using DifferentialEquations, BEFWM2, Distributions, ProgressMeter
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")

import Random.seed!

seed!(22)

A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0]
foodweb = FoodWeb(A)

# Parameter product
#
# Along correlation between environmental stochasticity and strength of
# stochasticity
ρ = -1.0:0.2:1.0
σₑ = 0.0:.2:1.5
rep = 1:10
Z_levels = [10:10:100;] # average predator-prey mass ratio
fr_types = [
    (h = 1.0, c = 0.0, type = "type II"),
    (h = 1.0, c = 1.0, type = "PI"),
    (h = 2.0, c = 0.0, type = "type III"),
]
names = (:ρ, :σₑ, :rep, :Z, :fr)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(ρ, σₑ, rep, Z_levels, fr_types))[:]

# Build the covariance matrix of the size of two times the number of species
S = BEFWM2.richness(foodweb)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc


println("Test with batch size = 5")
@elapsed testres = @showprogress pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 4000, dt = .1, corr_mat = vc), param[10000:10010], batch_size = 5)
println("Test with batch size = 1")
@elapsed testres2 =  @showprogress pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 4000, dt = .1, corr_mat = vc), param[10000:10010], batch_size = 1)

println("Real simulation")
@elapsed res = @showprogress pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 4000, dt = .1, corr_mat = vc), param, batch_size = 100)

df = DataFrame(res)

CSV.write("vasseur_fox__brose_res2.csv", df)
