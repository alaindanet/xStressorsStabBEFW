using Revise 
using BEFWM2
using Plots
using LinearAlgebra
using DifferentialEquations 
using DataFrames
using CSV

A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0] 
foodweb = FoodWeb(A)

# Functional response
## Preference of consumer 
myω = zeros(4, 4)
myω[:,1] = [0, 1, .98, 0]
myω[4,:] = [0, .92, 1 - .92, 0]
## Predator interference
myc = repeat([0], 4)
myc

# Biological rates
## Assimilitation efficiency or ingestion rates (in Vasseur & Fox, 2007)
## Efficiency is implicitly one for all consumers in Vasseur & Fox
mye = zeros(4, 4)
# Consumer
mye[:,1] = [0, 1, 1, 0] 
# Predator
mye[4,:] = [0, 1, 1, 0] 

# Metabolic rates and maximum ingestion rates (merged in J parameter in Vasseur
# & Fox)
J = [0, 0.8036, 0.7, .4]
## From McCann (1998):
x = [0, .40, .20, .08]
y = [0, 2.009, 3.50, 5.0]
x .* y .- J

bioener = BioenergeticResponse(foodweb,
                               h = 1,
                               # Half saturation-constant
                               B0 = [0, 0.16129, .9, .5],
                               # Consumer preference
                               ω = myω,
                               # Predator interference
                               c = myc
                              )

biorate = BioRates(foodweb,
        r = [1.0, 0, 0, 0],
        e = mye,
        x = x, 
        y = y
       )


params = ModelParameters(foodweb,
                functional_response = bioener,
                biorates = biorate,
                environment = Environment(foodweb, K = 1.0),
                env_stoch = EnvStoch(0.0)
               )

m = simulate(params, rand(4))

plot(m)


#
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

m = simulate(params, rand(4), diff_code_data=(mydBdt!, params))
plot(m)

function stochastic_process(dW, B, params, t)

    S = length(params.network.species)

    # Biomass dynamics have no stochasticity 
    for i in 1:2*S 
        dW[i] = 0.0
    end

    # Stochastic mortality variable for consumer is... Stochastic!!! 
    dW[S+2:S+3] .= params.env_stoch.σₑ 

end


# Stochastic strengh and correlation 
ρₑ = 0 

# Build the covariance matrix of the size of two times the number of species
S = BEFWM2.richness(foodweb)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc[S+2, S+3] = ρₑ
vc[S+3, S+2] = ρₑ
vc
 
# Generate the Wiener Process
W = CorrelatedWienerProcess(vc, 0.0, zeros(size(vc, 1)))

# define the starting values 
stoch_starting_val = [0; 0; 0; 0]
u0 = [rand(S); stoch_starting_val]
tspan = (0, 1000)

prob = SDEProblem(
           mydBdt!,
           stochastic_process,
           u0,
           tspan,
           params,
           noise = W
          )
sol = solve(prob;
      saveat = collect(0:1:1000),
      dt = .25,
      adaptive = false
      )
plot(sol, idxs = [1, 2, 3, 4])
plot(sol, idxs = 5:8)

foodweb_cv(sol, last = 100, idxs = [1, 2, 3, 4])

# Test with some noise 
params = ModelParameters(foodweb,
                functional_response = bioener,
                biorates = biorate,
                environment = Environment(foodweb, K = 1.0),
                env_stoch = EnvStoch(0.25)
               )
u0 = [rand(S); stoch_starting_val]
tspan = (0, 1000)
prob = SDEProblem(
           mydBdt!,
           stochastic_process,
           u0,
           tspan,
           params,
           noise = W
          )
sol = solve(prob;
      saveat = collect(0:1:1000),
      dt = .25,
      adaptive = false
      )
plot(sol, idxs = [1, 2, 3, 4])
plot(sol, idxs = 5:8)

foodweb_cv(sol, last = 100, idxs = [1, 2, 3, 4])
