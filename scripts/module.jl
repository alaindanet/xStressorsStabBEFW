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
@everywhere include("../src/get_modules.jl")

# include("src/stochastic_mortality_model.jl")
# include("src/sim.jl")
# include("src/interaction_strength.jl")

import Random.seed!

fw_module = get_fw_modules()

#sim_int_mat(A; σₑ = .5, alpha_ij = 0, Z = 100, h = 2.0, c = 1.0, K = 1.0, max = 50000, last = 25000, dt = 0.1, return_sol = false)
#
sim = map(x -> sim_int_mat(x;
			    ρ = 0, alpha_ij = 1.0,
			    σₑ = .5,
			    d = 0.1,
			    Z = 100, h = 2.0, c = 1.0, K = 3.0,
			    fun = stoch_d_dBdt!,
			    max = 1000, last = 100, dt = 0.1, return_sol = false),
	    fw_module)

rep = 1:20
fw_module_names = keys(fw_module)
ρ = [0.0, .5, 1.0]
alphaij = [0.0, .5, .8, 1]
sigma = [0, 0.5, .8]
Z = [5, 10, 100]
param_names = (:rep, :module, :rho, :alphaij,:sigma, :Z)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, fw_module_names, ρ, alphaij, sigma, Z))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, module_name = p.module),
                     sim_int_mat(fw_module[p.module];
                                 ρ = p.rho, alpha_ij = p.alphaij,
                                 σₑ = p.sigma,
                                 d = 0.1,
                                 Z = p.Z, h = 2.0, c = 1.0,
                                 r = 1.8, K = 20.0,
                                 fun = stoch_d_dBdt!,
                                 max = 5000, last = 1000, dt = 0.1, return_sol = false)
                    ),
                         param
                        )

df = DataFrame(sim)
Arrow.write("sim_module_dZ.arrow", df)
