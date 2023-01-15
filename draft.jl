using Revise
using BEFWM2
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
using JSON3
using CSV
include("src/minmax.jl")
include("src/interaction_strength.jl")
include("src/stochastic_mortality_model.jl")
include("src/sim.jl")

#simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
ti = simCS(0.1, 20, 100, 2.0, 1.0, 1.0, 5; max = 5000, last = 1000, dt = 0.1, return_sol = false)
plot(ti, idxs = collect(1:1:length(get_parameters(ti).network.species)))
cv(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))
biomass(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))
BEFWM2.filter_sim(ti)
get_parameters(ti)



# fw = FoodWeb([0 0 0; 1 0 0; 0 1 0], Z = 100)
fw = FoodWeb([0 0 0; 0 0 0; 0 0 0], Z = 100)
p = ModelParameters(fw,
                    functional_response = BioenergeticResponse(fw, h = 2, c = 1),
                    producer_competition = ProducerCompetition(fw; αij = .5),
                    env_stoch = EnvStoch(.5),
                    biorates = BioRates(fw; d = 0.05)
                   )
S = size(fw.A, 1)

stoch_starting_val = repeat([0], S)
u0 = [rand(S); stoch_starting_val]

# Make the stochastic matrix
corr_mat = zeros(S * 2, S * 2)
corr_mat .= 1.0
corr_mat[diagind(corr_mat)] .= 1.0
# Generate the Wiener Process
wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

prob = SDEProblem(
                  stoch_d_dBdt!,
                  gen_stochastic_process,
                  u0,
                  [0, 1000],
                  p,
                  noise = wiener
                 )
# Simulate
m = solve(prob;
      saveat = collect(0:1:1000),
      dt = .1,
      adaptive = false
     )
get_stab_fw(m; last = 100)
cv(m, last = 100, idxs = collect(1:1:S))
BEFWM2.synchrony(transpose(m[S+1:1:2*S, end-(100-1):end]))

biomass(m, last = 100, idxs = collect(1:1:S))
plot(m, idxs = collect(1:1:S))
plot(m, idxs = collect(S+1:1:2 * S))

sim_int_mat([0 0 0; 0 0 0; 1 1 0];
            ρ = 1.0, alpha_ij = 0,
            d = 0.1,
            σₑ = .5, Z = 100, h = 2.0, c = 1.0, K = 1.0,
            fun = stoch_d_dBdt!,
            max = 500, last = 100, dt = 0.1, return_sol = false)

