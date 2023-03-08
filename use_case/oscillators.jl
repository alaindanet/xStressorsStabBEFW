using EcologicalNetworksDynamics
using Plots
using LinearAlgebra
using DifferentialEquations
using DataFrames
using CSV




#########################################
#  Hajian-Forooshan (2020) 3 consumers  #
#########################################
include("../src/get_modules.jl")


# Define parameters
#
## Consumers
### Mortality rates for consumers
d_consumers = .1
### Handling time for consumers
ht = 3.0
### Attack rate
ar = .7
# Preference for secondary resources
B = .01

##
## Growth rate of producers
r = .3
## Carrying capacity
K = 1.0
## Competition among producers
α = .1

function get_parameters_kuramoto(;
        A = kuramoto_trophic_module()[1],
        comp_A = kuramoto_competition_module()[1],
        dcons = .1,
        ht = 3.0,
        ar = .7,
        B = .01,
        r = .3,
        K = 1.0,
        α = .1
    )

    # Preference matrix
    pref = [
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            1 B B 0 0 0;
            B 1 B 0 0 0;
            B B 1 0 0 0;
           ]
    # Set preference to 0 for non-interacting species pair
    mod_bool = Array{Bool}(A)
    pref[.!mod_bool] .= 0

    # Foodweb
    foodweb = FoodWeb(A)
    # Functional response
    classicresp = ClassicResponse(foodweb,
                                  # Hill exponent
                                  h = 1.0,
                                  # Consumer preference
                                  ω = pref,
                                  # Predator interference
                                  c = repeat([0], 6),
                                  # Handling time
                                  hₜ = ht,
                                  # Attack rates
                                  aᵣ = ar
                                 )
    # Biorates
    bio_rate = BioRates(foodweb,
                        d = [0, 0, 0, dcons, dcons, dcons],
                        r = [r, r, r, 0, 0, 0],
                        e = ones(6, 6),
                        x = 0.0,
                        y = 0.0
                       )
    # Environment
    env = Environment(foodweb, K = K)
    # Prod competition
    comp = [
            1 α α 0 0 0;
            α 1 α 0 0 0;
            α α 1 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
           ]
    comp_kuramoto_bool = Array{Bool}(comp_A)
    comp[.!comp_kuramoto_bool] .= 0
    prod_comp = ProducerCompetition(foodweb, α = comp)

    # ModelParameters
    ModelParameters(foodweb,
                        functional_response = classicresp,
                        biorates = bio_rate,
                        environment = env,
                        producer_competition = prod_comp
                       );
end

plot_storage = []
for i in 1:4

    p = get_parameters_kuramoto(
                                A = kuramoto_trophic_module()[i],
                                comp_A = kuramoto_competition_module()[i]
                               )
    m = simulate(p, [.05, .1, .35, .06, .11, .36],
                 tmax = 600,
                 callback = nothing
                );
    tmp_plot = plot(m, idxs = [4, 5, 6], legend = false)
    push!(plot_storage, tmp_plot)
end

myplot = plot(plot_storage[1], plot_storage[2], plot_storage[3], plot_storage[4])
png(myplot, "hajian_forooshan_vandermeer.png")

p = get_parameters_kuramoto(
                            A = kuramoto_trophic_module()[4],
                            comp_A = kuramoto_competition_module()[3],
                            B = .02, α = 0.1
                           )
m = simulate(p, [.05, .1, .35, .06, .11, .36],
             tmax = 2000,
             callback = nothing
            );
plot(m, idxs = [4, 5, 6])


plot_storage = []
for i in 1:4

    p = get_parameters_two_species_kuramoto()
    m = simulate(p, [.05, .1, .35, .06, .11, .36],
                 tmax = 600,
                 callback = nothing
                );
    tmp_plot = plot(m, idxs = [4, 5, 6], legend = false)
    push!(plot_storage, tmp_plot)
end

myplot = plot(plot_storage[1], plot_storage[2], plot_storage[3], plot_storage[4])


#######################
#  Vandermeer (2004)  #
#######################

A = [0 0 0 0; 0 0 0 0; 1 1 0 0; 1 1 0 0]
foodweb = FoodWeb(A)

# Functional response
## Preference of consumer
myω = zeros(4, 4)
# Consumer 1
B = 0.0015
myω[3,:] = [1, B, 0, 0]
myω[4,:] = [B, 1, 0, 0]

classicresp = ClassicResponse(foodweb,
                              h = 1.0,
                              # Consumer preference
                              ω = myω,
                              # Predator interference
                              c = repeat([0], 4),
                              # Handling time
                              hₜ = 1.3,
                              aᵣ = 2.0
                             )

bio_rate = BioRates(foodweb,
                    d = [0, 0, 0.1, 0.1],
                    r = [1.0, 1.0, 0, 0],
                    e = ones(4, 4),
                    x = 0.0,
                    y = 0.0
                   )

p = ModelParameters(foodweb,
                    functional_response = classicresp,
                    biorates = bio_rate,
                    environment = Environment(foodweb, K = 1.0),
                    producer_competition = ProducerCompetition(foodweb, αii = 0, αij = 0)
                   );

m = simulate(p, repeat([.5], 4));
plot(m)


##################
#  Bioenergetic  #
##################

#
foodweb2 = FoodWeb(A, Z = 10)
## Preference of consumer
myω = zeros(4, 4)
# Consumer 1
B = 0.10
#B = 0.00001
B = 0.0015
myω[3,:] = [1-B, B, 0, 0]
myω[4,:] = [B, 1-B, 0, 0]
p = ModelParameters(foodweb2,
                    functional_response = BioenergeticResponse(foodweb2, h = 1.0, ω = myω),
                    biorates = BioRates(foodweb2),
                    environment = Environment(foodweb2, K = 1.0),
                    producer_competition = ProducerCompetition(foodweb2, αii = 0.1, αij = 0.1)
                   );

# m = simulate(p, repeat([.5], 4));
m = simulate(p, repeat([.2,.3], 2));
m = simulate(p, rand(4), tmax = 1200);
plot(m, idxs = [3,4])
biomass(m)
