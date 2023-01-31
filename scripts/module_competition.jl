
##############################
#  Module and stochasticity  #
##############################


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

@everywhere import Pkg
@everywhere using DifferentialEquations, BEFWM2, Distributions, ProgressMeter, SparseArrays, LinearAlgebra, DataFrames, CSV, Arrow
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../src/interaction_strength.jl")

################################################################################
#                             Competition modules                              #
################################################################################

rep = 1:20
S = 1:10
ρ = [0.0, .5, 1.0]
alphaij = [0, .2, .5, .8, 1.0]
sigma = [.5, .9]
param_names = (:rep, :S, :rho, :alphaij,:sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, S, ρ, alphaij, sigma))[:]

sim_competition = @showprogress pmap(p -> merge(
                     (rep = p.rep, S = p.S),
                     sim_int_mat(zeros(p.S, p.S);
                                 ρ = p.rho, alpha_ij = p.alphaij,
                                 σₑ = p.sigma,
                                 d = 0.1,
                                 Z = 10, h = 2.0, c = 1.0,
                                 r = 1.8, K = 20.0,
                                 K_alpha_corrected = true,
                                 fun = stoch_d_dBdt!,
                                 max = 2000, last = 1000, dt = 0.1, return_sol = false)
                    ),
                                     param
                                    )

df_competition = DataFrame(sim_competition)
Arrow.write("sim_module_competition.arrow", df_competition)
