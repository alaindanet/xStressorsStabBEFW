# Test Thebault et al.
#

import Pkg
Pkg.activate("..")
using BEFWM2, LinearAlgebra, DifferentialEquations, DataFrames, CSV, Distributed, Distributions, ProgressMeter, BenchmarkTools
using Plots

include("src/stochastic_mortality_model.jl")
include("src/sim.jl")

import Random.seed!

seed!(22)

specialist = [
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
    ]
generalist = [
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     1 1 1 0 0 0;
     1 1 1 0 0 0;
     1 1 1 0 0 0;
             ]
foodweb = FoodWeb(specialist, Z = 100)
nprod = sum(foodweb.metabolic_class .== "producer")
param = ModelParameters(foodweb,
                       environment = Environment(foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(foodweb, h = 1)
                      )
spe_sim = simulate(param, rand(length(foodweb.species)), tmax = 500, callback = nothing)
foodweb_cv(spe_sim, last = 100)
plot(spe_sim)


foodweb = FoodWeb(generalist, Z = 100)
nprod = sum(foodweb.metabolic_class .== "producer")
param = ModelParameters(foodweb,
                       environment = Environment(foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(foodweb, h = 1)
                      )
gen_sim = simulate(param, rand(length(foodweb.species)), tmax = 500, callback = nothing)
foodweb_cv(gen_sim, last = 100)
plot(gen_sim)


foodweb_cv(spe_sim, last = 50)
foodweb_cv(gen_sim, last = 50)

# Add environmental stochasticity
#
#
function stochastic_process(dW, B, params, t)

    S = length(params.network.species)

    # Biomass dynamics have no stochasticity
    for i in 1:2*S
        dW[i] = 0.0
    end

    # Stochastic mortality variable for consumer is... Stochastic!!!
    dW[S+1:2*S] .= params.env_stoch.σₑ

end

#
S = BEFWM2.richness(specialist)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc
stoch_starting_val = [0; 0; 0; 0; 0; 0]
wiener = CorrelatedWienerProcess(vc, 0.0, zeros(size(vc, 1)))

spe_foodweb = FoodWeb(specialist, Z = 100)
nprod = sum(spe_foodweb.metabolic_class .== "producer")
param = ModelParameters(spe_foodweb,
                       environment = Environment(spe_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(spe_foodweb, h = 2, c = 1.0),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# Generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
spe_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(spe_sim_stoch, last = 50, idxs = [1:6;])
plot(spe_sim_stoch, idxs = [1:6;])


gen_foodweb = FoodWeb(generalist, Z = 100)
nprod = sum(gen_foodweb.metabolic_class .== "producer")
param = ModelParameters(gen_foodweb,
                       environment = Environment(gen_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(gen_foodweb, h = 2, c = 1.0),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# Generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
gen_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(gen_sim_stoch, last = 50, idxs = [1:6;])
plot(gen_sim_stoch, idxs = [1:6;])

foodweb_cv(gen_sim_stoch, last = 2500, idxs = [1:6;])
foodweb_cv(spe_sim_stoch, last = 2500, idxs = [1:6;])

foodweb_cv(gen_sim_stoch, last = 2500, idxs = [4:6;])
foodweb_cv(spe_sim_stoch, last = 2500, idxs = [4:6;])


# Same but with a top generalist consumer
#
top_specialist = [
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     1 0 0 0 0 0 0;
     0 1 0 0 0 0 0;
     0 0 1 0 0 0 0;
     0 0 0 1 1 1 0;
    ]
top_generalist = [
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     1 1 1 0 0 0 0;
     1 1 1 0 0 0 0;
     1 1 1 0 0 0 0;
     0 0 0 1 1 1 0;
             ]

S = BEFWM2.richness(top_specialist)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc
stoch_starting_val = zeros(S)
wiener = CorrelatedWienerProcess(vc, 0.0, zeros(size(vc, 1)))

top_spe_foodweb = FoodWeb(top_specialist, Z = 10)
nprod = sum(top_spe_foodweb.metabolic_class .== "producer")
param = ModelParameters(top_spe_foodweb,
                       environment = Environment(top_spe_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(top_spe_foodweb, h = 2, c = 1),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# Generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
top_spe_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(top_spe_sim_stoch, last = 50, idxs = [1:6;])
plot(top_spe_sim_stoch, idxs = [1:6;])


top_gen_foodweb = FoodWeb(top_generalist, Z = 10)
nprod = sum(top_gen_foodweb.metabolic_class .== "producer")
param = ModelParameters(top_gen_foodweb,
                       environment = Environment(top_gen_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(top_gen_foodweb, h = 2, c = 1),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# top_generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
top_gen_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(top_gen_sim_stoch, last = 2500, idxs = [1:6;])
plot(top_gen_sim_stoch, idxs = [1:6;])

foodweb_cv(top_gen_sim_stoch, last = 2500, idxs = [1:6;])
foodweb_cv(top_spe_sim_stoch, last = 2500, idxs = [1:6;])

foodweb_cv(top_gen_sim_stoch, last = 2500, idxs = [4:6;])
foodweb_cv(top_spe_sim_stoch, last = 2500, idxs = [4:6;])


# Same but with a top generalist consumer
#
top_omn_specialist = [
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     1 0 0 0 0 0 0;
     0 1 0 0 0 0 0;
     0 0 1 0 0 0 0;
     1 1 1 1 1 1 0;
    ]
top_omn_generalist = [
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     1 1 1 0 0 0 0;
     1 1 1 0 0 0 0;
     1 1 1 0 0 0 0;
     1 1 1 1 1 1 0;
             ]

S = BEFWM2.richness(top_omn_specialist)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc
stoch_starting_val = zeros(S)
wiener = CorrelatedWienerProcess(vc, 0.0, zeros(size(vc, 1)))

top_omn_spe_foodweb = FoodWeb(top_omn_specialist, Z = 10)
nprod = sum(top_omn_spe_foodweb.metabolic_class .== "producer")
param = ModelParameters(top_omn_spe_foodweb,
                       environment = Environment(top_omn_spe_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(top_omn_spe_foodweb, h = 2, c = 1),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# Generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
top_omn_spe_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(top_omn_spe_sim_stoch, last = 50, idxs = [1:6;])
plot(top_omn_spe_sim_stoch, idxs = [1:6;])


top_omn_gen_foodweb = FoodWeb(top_omn_generalist, Z = 10)
nprod = sum(top_omn_gen_foodweb.metabolic_class .== "producer")
param = ModelParameters(top_omn_gen_foodweb,
                       environment = Environment(top_omn_gen_foodweb,K = 3/nprod),
                       functional_response = BioenergeticResponse(top_omn_gen_foodweb, h = 2, c = 1),
                       env_stoch = EnvStoch(.5)
                      )
u0 = [rand(S); stoch_starting_val]

# top_omn_generate the Wiener Process

prob = SDEProblem(
                  mydBdt!,
                  stochastic_process,
                  u0,
                  [0, 5000],
                  param,
                  noise = wiener
                 )
# Simulate
top_omn_gen_sim_stoch = solve(prob;
      saveat = collect(0:1:5000),
      dt = .1,
      adaptive = false
     )
foodweb_cv(top_omn_gen_sim_stoch, last = 2500, idxs = [1:6;])
plot(top_omn_gen_sim_stoch, idxs = [1:6;])

foodweb_cv(top_omn_gen_sim_stoch, last = 2500, idxs = [1:7;])
foodweb_cv(top_omn_spe_sim_stoch, last = 2500, idxs = [1:7;])

foodweb_cv(top_omn_gen_sim_stoch, last = 2500, idxs = [4:6;])
foodweb_cv(top_omn_spe_sim_stoch, last = 2500, idxs = [4:6;])


foodweb_cv(top_omn_gen_sim_stoch, last = 2500, idxs = [1:7;])
foodweb_cv(top_gen_sim_stoch, last = 2500, idxs = [1:7;])

foodweb_cv(top_omn_spe_sim_stoch, last = 2500, idxs = [1:7;])
foodweb_cv(top_spe_sim_stoch, last = 2500, idxs = [1:7;])


function sim_module(A; h = 2, c = 1, σ = .5, z = 10, k_global = 3, tmax = 5000, last = 2500, return_sol = false)

    fw = FoodWeb(A, Z = z)

    S = BEFWM2.richness(fw)
    vc = zeros(S * 2, S * 2)
    # Only the consumer have a variance
    # vc[diagind(vc)][S+2:S+3] .= σₑ
    vc[diagind(vc)] .= 1.0
    vc

    stoch_starting_val = zeros(S)
    wiener = CorrelatedWienerProcess(vc, 0.0, zeros(size(vc, 1)))
    nprod = sum(fw.metabolic_class .== "producer")
    param = ModelParameters(fw,
                            environment = Environment(fw,K = k_global/nprod),
                            functional_response = BioenergeticResponse(fw, h = h, c = c),
                            env_stoch = EnvStoch(σ)
                           )
    u0 = [rand(S); stoch_starting_val]

    # Generate the Wiener Process

    prob = SDEProblem(
                      mydBdt!,
                      stochastic_process,
                      u0,
                      [0, tmax],
                      param,
                      noise = wiener
                     )
    # Simulate
    sol = solve(prob;
          saveat = collect(0:1:tmax),
          dt = .1,
          adaptive = false
         )
    if return_sol
        return sol
    end
    foodweb_cv(sol, last = last, idxs = [1:S;])
end


# Parameter product
#
#
rep = 1:10
foodwebs = [
            (name = "spe", A = specialist),
            (name = "gen", A = generalist),
            (name = "top_spe", A = top_specialist),
            (name = "top_gen", A = top_generalist),
            (name = "top_omn_spe", A = top_omn_specialist),
            (name = "top_omn_gen", A = top_omn_generalist),

           ]
names = (:rep, :foodweb)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(rep, foodwebs))[:]

sim = map(p -> merge(
                     (rep = p.rep, fw = p.foodweb.name),
                     sim_module(p.foodweb.A; h = 2, c = 1, σ = .5, z = 10, k_global = 3, tmax = 5000, last = 2500, return_sol = false)
                    ),
          param)

df = DataFrame(sim)
CSV.write("generalism_stoch.csv", df)
