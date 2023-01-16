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
@everywhere include("src/stochastic_mortality_model.jl")
@everywhere include("src/sim.jl")
@everywhere include("src/interaction_strength.jl")

# include("src/stochastic_mortality_model.jl")
# include("src/sim.jl")
# include("src/interaction_strength.jl")

import Random.seed!

fw_module = (
          # prod = [0 0; 0 0],
          chain1 = [0 0; 1 0],
          chain2 = [0 0 0; 1 0 0; 0 1 0],
          chain3 = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0],
          chain4 = [0 0 0 0;
                    1 0 0 0;
                    0 1 0 0;
                    0 0 1 0],
          cons = [0 0 0; 0 0 0; 1 1 0],
          cons2_spe = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0],
          cons2_spe_top = [0 0 0 0 0;
                           0 0 0 0 0;
                           1 0 0 0 0;
                           0 1 0 0 0;
                           0 0 1 1 0],
          cons3_spe = [0 0 0 0 0 0;
                       0 0 0 0 0 0;
                       1 0 0 0 0 0;
                       0 1 0 0 0 0;
                       0 0 1 0 0 0;
                       0 0 0 1 0 0],
          cons2_gen = [0 0 0 0; 0 0 0 0; 1 1 0 0; 1 1 0 0],
          cons2_gen_top = [0 0 0 0 0;
                           0 0 0 0 0;
                           1 1 0 0 0;
                           1 1 0 0 0;
                           0 0 1 1 0],
          cons3_gen = [0 0 0 0 0 0;
                       0 0 0 0 0 0;
                       1 1 0 0 0 0;
                       1 1 0 0 0 0;
                       0 0 1 1 0 0;
                       0 0 1 1 0 0],
         )

#sim_int_mat(A; σₑ = .5, alpha_ij = 0, Z = 100, h = 2.0, c = 1.0, K = 1.0, max = 50000, last = 25000, dt = 0.1, return_sol = false)
#
sim = map(x -> sim_int_mat(x;
                                 ρ = 0, alpha_ij = 1.0,
                                 σₑ = .5,
                                 d = 0.1,
                                 Z = 100, h = 2.0, c = 1.0, K = 3.0,
                                 fun = stoch_d_dBdt!,
                                 max = 1000, last = 100, dt = 0.1, return_sol = false)
                     , fw_module)

rep = 1:20
fw_module_names = keys(fw_module)
ρ = [0.0, .5, 1.0]
alphaij = [0.0, .5, 1]
sigma = [0.5, 0.9]
param_names = (:rep, :module, :rho, :alphaij,:sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, fw_module_names, ρ, alphaij, sigma))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, module_name = p.module),
                     sim_int_mat(fw_module[p.module];
                                 ρ = p.rho, alpha_ij = p.alphaij,
                                 σₑ = p.sigma,
                                 d = 0.1,
                                 Z = 10, h = 2.0, c = 1.0,
                                 r = 1.8, K = 20.0,
                                 fun = stoch_d_dBdt!,
                                 max = 5000, last = 1000, dt = 0.1, return_sol = false)
                    ),
                         param
                        )

df = DataFrame(sim)
Arrow.write("res/sim_module_d.arrow", df)

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
Arrow.write("res/sim_module_competition.arrow", df_competition)