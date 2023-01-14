##############################
#  Module and stochasticity  #
##############################


import Pkg
Pkg.activate(".")
using BEFWM2, LinearAlgebra, DifferentialEquations, DataFrames, CSV, Distributed, Distributions, ProgressMeter, BenchmarkTools
using Plots, SparseArrays
using Arrow

include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/interaction_strength.jl")

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

map(x -> FoodWeb(x), fw_module)

foodweb = FoodWeb(specialist, Z = 100)
nprod = sum(foodweb.metabolic_class .== "producer")
param = ModelParameters(foodweb,
                       environment = Environment(foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(foodweb, h = 1)
                      )

#sim_int_mat(A; σₑ = .5, alpha_ij = 0, Z = 100, h = 2.0, c = 1.0, K = 1.0, max = 50000, last = 25000, dt = 0.1, return_sol = false)
sim = map(x -> sim_int_mat(x, .5, .5; Z = 100, h = 2.0, c = 1.0, K = 3, max = 5000, last = 1000, dt = 0.1, return_sol = false), fw_module)

rep = 1:20
fw_module_names = keys(fw_module)
Z = [5, 10, 100]
sigma = 0.5:1
param_names = (:rep, :module, :Z, :sigma)
param = map(p -> (;Dict(k => v for (k, v) in zip(param_names, p))...), Iterators.product(rep, fw_module_names, alpha_ij, sigma))[:]

sim = @showprogress pmap(p -> merge(
                     (rep = p.rep, alphaij = p.alphaij),
                     sim_int_mat(fw_module[p.module];
                                 σₑ = p.sigma, alpha_ij = 1.0,
                                 Z = p.Z, h = 2.0, c = 1.0, K = 3,
                                 max = 5000, last = 1000, dt = 0.1, return_sol = false
                                )
                    ),
           param
          )


df = DataFrame(sim)
Arrow.write("res/sim_moduleZ.arrow", df)

# Some playing with Arrow
# table2 = Arrow.Table("test_arrow.arrow")
# dt = DataFrame(table2)
# df.int_strength[240]
# dt.int_strength[240]
# reshape(dt.int_strength[240], 5, 5)

# Arrows
