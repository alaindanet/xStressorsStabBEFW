using Revise
using EcologicalNetworksDynamics
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Distributions
using Statistics
using DataFrames
using Plots
using Debugger
using CSV
using Arrow
include("src/minmax.jl")
include("src/interaction_strength.jl")
#include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/plot.jl")
include("src/get_modules.jl")

const fractional_digits = 4

foodweb = FoodWeb([0 0 0; 1 0 0; 0 1 0]);
params = ModelParameters(foodweb,
                         biorates = BioRates(foodweb, d = .1),
                         env_stoch = EnvStoch(.5)
                        );
B0 = [0.5, 0.5, 0.5];
sol = simulate(params, B0;
               rho = 1,
               dt = .1,
               tmax = 5000,
               extinction_threshold = 1e-5,
               verbose = false
              );
mat_rounded = round.(sol[:,:], digits = 5)
mat = sol[:,:]
coefficient_of_variation(sol)
cv_sp(sol)
ti = empirical_interaction_strength(sol, params, last = 100)
trophic_structure(sol, last = 100)

richness(sol)
species_persistence(sol)
ti = rand(2, 2, 3)
ti.round()
round.(ti; digits = 5)
ti[:,:, 1]
richness(params.network)

ti = simCS(.1, 20, ρ = 1.0, max = 200, last = 100)
ti.time_species
ti.time_stoch
varinfo(r"ti")



#########
#  Sim  #
#########

#simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
using Distributions



ti = simCS(0.1, 20;
           d = .1,
           Z = 100, h = 2.0, c = 0.0,
           σₑ = 1.0, K = 5, max = 5000,
           last = 100, dt = 0.1, return_sol = true
          )

plot(ti, idxs = collect(1:1:length(get_parameters(ti).network.species)))
cv(ti, last = 100, idxs = collect(1:1:length(get_parameters(ti).network.species)))
biomass(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))



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
EcologicalNetworksDynamics.synchrony(transpose(m[S+1:1:2*S, end-(100-1):end]))

biomass(m, last = 100, idxs = collect(1:1:S))
plot(m, idxs = collect(1:1:S))
plot(m, idxs = collect(S+1:1:2 * S))

sim_int_mat([0 0 0; 0 0 0; 1 1 0];
            ρ = 1.0, alpha_ij = 0,
            d = 0.1,
            σₑ = .5, Z = 100, h = 2.0, c = 1.0, K = 1.0,
            fun = stoch_d_dBdt!,
            max = 500, last = 100, dt = 0.1, return_sol = false)



###########
#  Motif  #
###########

fw = FoodWeb(nichemodel, 20, C = .05)
fw = FoodWeb(nichemodel, 80, C = .05)
EcologicalNetworks.find_motif(fw.A, unipartitemotifs().S1)
mot = find_motif(UnipartiteNetwork(fw.A), unipartitemotifs().S1) |> length
map(x -> find_motif(UnipartiteNetwork(fw.A), x) |> length, unipartitemotifs())

# Diamond
map(x -> find_motif(UnipartiteNetwork([0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0] .> 0), x) |> length, unipartitemotifs())

##########
#  Plot  #
##########

webplot(get_fw_modules()[10]; consasrow = true)
webplot([0 0 0; 1 1 0; 1 1 0]; consasrow = true)
ti = map(x -> webplot(x; consasrow = true), get_fw_modules())
length(ti)
plot(ti..., layout = (4, 4))
plot(ti[1], ti[2], ti[3], layout = (1, 3))
plot((ti[i] for i in 1:length(ti))..., layout = (4, 3))

