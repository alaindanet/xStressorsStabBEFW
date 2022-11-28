function mydBdt!(dB, B, params::ModelParameters, t)

    # Set up - Unpack parameters
    S = length(params.network.species)
    response_matrix = params.functional_response(B[1:S], params.network)
    r = params.biorates.r # vector of intrinsic growth rates
    K = params.environment.K # vector of carrying capacities
    network = params.network
    θ = 1

    # Compute ODE terms for each species
    for i in 1:S
        growth = BEFWM2.logisticgrowth(i, B[1:S], r[i], K[i], network)
        eating, being_eaten = BEFWM2.consumption(i, B[1:S], params, response_matrix)
        # Metabolic loss as basal mortality times exponentional stochastic noise
        metabolism_loss = params.biorates.x[i] * exp(B[i+S]) * B[i]
        net_growth_rate = growth + eating - metabolism_loss
        # net_growth_rate = BEFWM2.effect_competition(net_growth_rate, i, B, network)
        dB[i] = net_growth_rate - being_eaten
    end

    # Compute stochastic mortality variable for all species
    for i in S+1:2*S
        # Ornstein-Uhlenbeck process
        dB[i] = θ * (0 - B[i])
    end
end

function stochastic_process(dW, B, params, t)

    S = length(params.network.species)

    # Biomass dynamics have no stochasticity
    for i in 1:2*S
        dW[i] = 0.0
    end

    # Stochastic mortality variable for consumer is... Stochastic!!!
    dW[S+2:S+3] .= params.env_stoch.σₑ

end
