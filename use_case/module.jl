##############################
#  Module and stochasticity  #
##############################


import Pkg
using DifferentialEquations, EcologicalNetworksDynamics, ProgressMeter, SparseArrays, LinearAlgebra, DataFrames, CSV, Arrow
using Distributions
include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/interaction_strength.jl")

# include("src/stochastic_mortality_model.jl")
# include("src/sim.jl")
# include("src/interaction_strength.jl")

import Random.seed!

fw_module = (
          prod = [0 0; 0 0],
          chain1 = [0 0; 1 0],
          chain2 = [0 0 0; 1 0 0; 0 1 0],
          chain3 = [0 0 0 0;
                    1 0 0 0;
                    0 1 0 0;
                    0 0 1 0],
          chain4 = [0 0 0 0 0;
                    1 0 0 0 0;
                    0 1 0 0 0;
                    0 0 1 0 0;
                    0 0 0 1 0],
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

module_fw = map(x -> FoodWeb(x), fw_module)

foodweb = FoodWeb(module_fw.cons, Z = 100)
nprod = sum(foodweb.metabolic_class .== "producer")
param = ModelParameters(foodweb,
                       environment = Environment(foodweb, K = 3/nprod),
                       functional_response = BioenergeticResponse(foodweb, h = 1)
                      )

#sim_int_mat(A; σₑ = .5, alpha_ij = 0, Z = 100, h = 2.0, c = 1.0, K = 1.0, max = 50000, last = 25000, dt = 0.1, return_sol = false)
#
#
alphaij = 0.0
sim = map(x -> sim_int_mat(x;
                           ρ = 0, alpha_ij = alphaij,
                           σₑ = 0.1,
                           d = 0.05,
                           Z = 100, h = 2.0, c = 1.0,
                           r = 1.0, K = 10.0,
                           K_alpha_corrected = true,
                           dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                           max = 2000, last = 100, dt = 0.1,
                           return_sol = false),
          fw_module);

df = DataFrame(sim)
df[:, :name] = collect(keys(fw_module))

df[:, [:name,:richness, :stab_com, :avg_cv_sp, :sync]]

rep = 1:10
fw_module_names = keys(fw_module)
ρ = [0.0, .5, 1.0]
alphaij = [0, .2, .5, .8, 1.0]
sigma = 0.5:.5:1.0
param_names = (:rep, :module, :rho, :alphaij,:sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, fw_module_names, ρ, alphaij, sigma))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, module_name = p.module),
                     sim_int_mat(fw_module[p.module];
                                 ρ = p.rho, alpha_ij = p.alphaij,
                                 σₑ = p.sigma,
                                 d = 0.1,
                                 Z = 100, h = 2.0, c = 1.0, K = 3.0,
                                 dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                                 max = 1000, last = 100, dt = 0.1, return_sol = false)
                    ),
                         param
                        )

df = DataFrame(sim)
Arrow.write("../res/sim_module_d.arrow", df)

FoodWeb(zeros(2,2), quiet = true)


rep = 1:10
fw_module_names = keys(fw_module)
ρ = [0.0, .5, 1.0]
alphaij = [0, .2, .5, .8, 1.0]
sigma = 0.5:.5:1.0
param_names = (:rep, :module, :rho, :alphaij,:sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, fw_module_names, ρ, alphaij, sigma))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, module_name = p.module),
                     sim_int_mat(fw_module[p.module];
                                 ρ = p.rho, alpha_ij = p.alphaij,
                                 σₑ = p.sigma,
                                 d = 0.1,
                                 Z = 100, h = 2.0, c = 1.0, K = 3.0,
                                 fun = stoch_d_dBdt!,
                                 max = 1000, last = 100, dt = 0.1, return_sol = false)
                    ),
                         param
                        )
sim_competition = map(x -> sim_int_mat(zeros(x, x);
                           ρ = 0, alpha_ij = .2,
                           σₑ = .1,
                           d = 0.1,
                           Z = 10, h = 2.0, c = 1.0,
                           r = 1.0, K = 10.0,
                           K_alpha_corrected = true,
                           dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                           max = 2000, last = 100, dt = 0.1,
                           return_sol = false),
          1:15);

sim_compet_df = DataFrame(sim_competition)
sim_compet_df[:, [:richness, :stab_com, :avg_cv_sp, :sync]]

df = DataFrame(sim)
