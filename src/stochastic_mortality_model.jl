function mydBdt!(dB, B, params::ModelParameters, t)

    # params, extinct_sp = p # unpack input

    # Set up - Unpack parameters
    S = richness(params.network)
    response_matrix = params.functional_response(B[1:S], params.network)
    r = params.biorates.r # vector of intrinsic growth rates
    d = params.biorates.d # vector of natural death
    K = params.environment.K # vector of carrying capacities
    α = params.producer_competition.α # matrix of producer competition
    network = params.network
    θ = 1


    for i in 1:S
        # sum(α[i, :] .* B)) measures competitive effects (s)
        growth = EcologicalNetworksDynamics.logisticgrowth(i, B, r[i], K[i], sum(α[i, :] .* B[1:S]), network)
        eating, being_eaten = EcologicalNetworksDynamics.consumption(i, B, params, response_matrix)
        metabolism_loss = params.biorates.x[i] * exp(B[i+S]) * B[i]
        natural_death = EcologicalNetworksDynamics.natural_death_loss(i, B, params)
        net_growth_rate = growth + eating - metabolism_loss
        # net_growth_rate = effect_competition(net_growth_rate, i, B, network)
        dB[i] = net_growth_rate - being_eaten - natural_death
    end

    # Compute stochastic mortality variable for all species
    for i in S+1:2*S
        # Ornstein-Uhlenbeck process
        dB[i] = θ * (0 - B[i])
    end

    # Avoid zombie species by forcing extinct biomasses to zero.
    # https://github.com/BecksLab/EcologicalNetworksDynamics/issues/65
    # for sp in keys(extinct_sp)
        # B[sp] = 0.0
    # end
end

function stoch_d_dBdt!(dB, B, params::ModelParameters, t)

    # params, extinct_sp = p # unpack input

    # Set up - Unpack parameters
    S = richness(params.network)
    response_matrix = params.functional_response(B[1:S], params.network)
    r = params.biorates.r # vector of intrinsic growth rates
    d = params.biorates.d # vector of natural death
    K = params.environment.K # vector of carrying capacities
    α = params.producer_competition.α # matrix of producer competition
    network = params.network
    θ = 1

    # Compute ODE terms for each species
    for i in 1:S
        # sum(α[i, :] .* B)) measures competitive effects (s)
        growth = EcologicalNetworksDynamics.logisticgrowth(i, B, r[i], K[i], sum(α[i, :] .* B[1:S]), network)
        eating, being_eaten = EcologicalNetworksDynamics.consumption(i, B, params, response_matrix)
        metabolism_loss = EcologicalNetworksDynamics.metabolic_loss(i, B, params)
        natural_death = d[i] * exp(B[i+S]) * B[i] # Stochastic natural death
        net_growth_rate = growth + eating - metabolism_loss
        # net_growth_rate = effect_competition(net_growth_rate, i, B, network)
        dB[i] = net_growth_rate - being_eaten - natural_death
    end

    # Avoid zombie species by forcing extinct biomasses to zero.
    # https://github.com/BecksLab/EcologicalNetworksDynamics/issues/65
    # for sp in keys(extinct_sp)
        # B[sp] = 0.0
    # end
    #
    # Compute stochastic mortality variable for all species
    for i in S+1:2*S
        # Ornstein-Uhlenbeck process
        dB[i] = θ * (0 - B[i])
    end

end

function stoch_m_dBdt!(dB, B, params::ModelParameters, t)

    # params, extinct_sp = p # unpack input

    # Set up - Unpack parameters
    S = richness(params.network)
    response_matrix = params.functional_response(B[1:S], params.network)
    r = params.biorates.r # vector of intrinsic growth rates
    d = params.biorates.d # vector of natural death
    K = params.environment.K # vector of carrying capacities
    α = params.producer_competition.α # matrix of producer competition
    network = params.network
    θ = 1

    # Compute ODE terms for each species
    for i in 1:S
        # sum(α[i, :] .* B)) measures competitive effects (s)
        growth = EcologicalNetworksDynamics.logisticgrowth(i, B, r[i], K[i], sum(α[i, :] .* B[1:S]), network)
        eating, being_eaten = EcologicalNetworksDynamics.consumption(i, B, params, response_matrix)
        metabolism_loss = EcologicalNetworksDynamics.metabolic_loss(i, B, params) * exp(B[i+S])
        natural_death = d[i] * exp(B[i+S]) * B[i] # Stochastic natural death
        net_growth_rate = growth + eating - metabolism_loss
        # net_growth_rate = effect_competition(net_growth_rate, i, B, network)
        dB[i] = net_growth_rate - being_eaten - natural_death
    end

    # Avoid zombie species by forcing extinct biomasses to zero.
    # https://github.com/BecksLab/EcologicalNetworksDynamics/issues/65
    # for sp in keys(extinct_sp)
        # B[sp] = 0.0
    # end
    #
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

function gen_stochastic_process(dW, B, params, t)

    S = length(params.network.species)

    # Biomass dynamics have no stochasticity
    dW .= 0.0

    # Stochastic mortality variable for consumer is... Stochastic!!!
    dW[S+1:2*S] .= params.env_stoch.σₑ

end
