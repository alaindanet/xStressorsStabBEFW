##############################
#  Module and stochasticity  #
##############################


import Pkg
Pkg.activate("..")
using BEFWM2, LinearAlgebra, DifferentialEquations, DataFrames, CSV, Distributed, Distributions, ProgressMeter, BenchmarkTools
using Plots

include("src/stochastic_mortality_model.jl")
include("src/sim.jl")

import Random.seed!

fw_module = (
          prod = [0 0; 0 0],
          chain1 = [0 0; 1 0],
          cons = [0 0 0; 0 0 0; 1 1 0],
          cons_spe = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0],
          cons_gen = [0 0 0 0; 0 0 0 0; 1 1 0 0; 1 1 0 0],
          cons_spe_top = [0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 1 0],
          cons_gen_top = [0 0 0 0 0; 0 0 0 0 0; 1 1 0 0 0; 1 1 0 0 0; 0 0 1 1 0],
         )

foodweb = FoodWeb(specialist, Z = 100)
nprod = sum(foodweb.metabolic_class .== "producer")
param = ModelParameters(foodweb,
                       environment = Environment(foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(foodweb, h = 1)
                      )

function simCS(A, σₑ, alpha_ij; Z, h, c, K, max = 50000, last = 25000, dt = 0.1, return_sol = false)

    fw = FoodWeb(A, Z = Z)

    # Parameters of the functional response
    funcrep = BioenergeticResponse(fw; h = h, c = c)
    # generate the model parameters
    p = ModelParameters(fw;
                        functional_response = funcrep,
                        environment = Environment(fw, K = K),
                        producer_competition = ProducerCompetition(fw, αii = 1.0, αij = alpha_ij),
                        env_stoch = EnvStoch(σₑ)
                       )
    stoch_starting_val = repeat([0], S)
    u0 = [rand(S); stoch_starting_val]

    # Make the stochastic matrix
    corr_mat = Diagonal(repeat([1.0], S * 2))
    # Generate the Wiener Process
    wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

    prob = SDEProblem(
                      mydBdt!,
                      gen_stochastic_process,
                      u0,
                      [0, max],
                      p,
                      noise = wiener
                     )
    # Simulate
    timing = @elapsed m = try
        solve(prob;
              saveat = collect(0:1:max),
              dt = dt,
              adaptive = false
             )

    catch
        (t = 0, x = missing)
    end

    if return_sol
        return m
    end

    if length(m.t) == max + 1
        cv = foodweb_cv(m, last = last, idxs = collect(1:1:S))
        bm = total_biomass(m, last = last, idxs = collect(1:1:S))
        int_strength = empirical_interaction_strength(m, p, last = last)
        non_zero_int = [i for i in int_strength if i != 0]
        max_int = max_interaction_strength(p)
        non_zero_max_int = [i for i in max_int if i != 0]
        tlvl = trophic_levels(fw)

        out = (
               richness          = S,
               ct                = C,
               Z                 = Z,
               env_stoch         = σₑ,
               sim_timing        = timing,
               stab_com          = 1 / cv.cv_com,
               avg_cv_sp         = cv.avg_cv_sp,
               sync              = cv.synchrony,
               total_biomass     = bm,
               bm_sp             = cv.bm_sp,
               cv_sp             = cv.cv_sp,
               tlvl              = tlvl,
               w_avg_tlvl        = sum(tlvl .* (cv.bm_sp ./ sum(cv.bm_sp))),
               max_tlvl          = maximum(tlvl),
               int_strength      = int_strength,
               avg_int_strength  = mean(non_zero_int),
               max_int_strength  = maximum(non_zero_int),
               min_int_strength  = minimum(non_zero_int),
               gini_int_strength = gini(non_zero_int),
               max_int           = max_int,
               avg_max_int       = mean(non_zero_max_int),
               max_max_int       = maximum(non_zero_max_int),
               min_max_int       = minimum(non_zero_max_int),
               gini_max_int      = gini(non_zero_max_int)
              )
    else
        out = (
               richness          = S,
               ct                = C,
               Z                 = Z,
               env_stoch         = σₑ,
               sim_timing        = timing,
               stab_com          = missing,
               avg_cv_sp         = missing,
               sync              = missing,
               total_biomass     = missing,
               bm_sp             = missing,
               cv_sp             = missing,
               tlvl              = missing,
               w_avg_tlvl        = missing,
               max_tlvl          = missing,
               int_strength      = missing,
               avg_int_strength  = missing,
               max_int_strength  = missing,
               min_int_strength  = missing,
               gini_int_strength = missing,
               max_int           = missing,
               avg_max_int       = missing,
               max_max_int       = missing,
               min_max_int       = missing,
               gini_max_int      = missing
              )
    end
    out
end
