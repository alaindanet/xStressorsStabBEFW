using Revise
using BEFWM2
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
include("src/minmax.jl")

Iₖ = ωₖ * x * y / B₀

######
# Define your own ModelParameters with an additional mortality term 
#####

# Define the composite type CritDeath
mutable struct CritDeath 
    d::Float64
end

# Add a field crit_death to ModelParameters by defining a new composite type of subtype Params
mutable struct CustomModelParameters <: Params
    network::EcologicalNetwork
    biorates::BioRates
    environment::Environment
    functional_response::FunctionalResponse
    crit_death::CritDeath
end


# Define the function of my custom ModelParameters
function CustomModelParameters(
    network::EcologicalNetwork;
    biorates::BioRates = BioRates(network),
    environment::Environment = Environment(network),
    functional_response::FunctionalResponse = BioenergeticResponse(network),
    crit_death::CritDeath = CritDeath(0.2)
)
    CustomModelParameters(network, biorates, environment, functional_response, crit_death)
end

# My custom dBdt! contains B[i] * params.crit_death.d in metabolic losses
# It takes a object of type Params in input but consumption() should also take a Params input instead of  #ModelParameters type of input
function CustomdBdt!(dB, B, params::Params, t)

    # Set up - Unpack parameters
    S = richness(params.network)
    response_matrix = params.functional_response(B, params.network)
    r = params.biorates.r # vector of intrinsic growth rates
    K = params.environment.K # vector of carrying capacities
    network = params.network

    # Compute ODE terms for each species
    for i in 1:S
        growth = logisticgrowth(i, B, r[i], K[i], network)
        eating, being_eaten = consumption(i, B, params, response_matrix)
        metabolism_loss = metabolic_loss(i, B, params) + B[i] * params.crit_death.d
        net_growth_rate = growth + eating - metabolism_loss
        net_growth_rate = effect_competition(net_growth_rate, i, B, network)
        dB[i] = net_growth_rate - being_eaten
    end
end

# Do the simulation
A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0] 
foodweb = FoodWeb(A)

params = CustomModelParameters(foodweb, crit_death = CritDeath(0.2))
m = simulate(params, rand(4), diff_code_data=(CustomdBdt!, params))

