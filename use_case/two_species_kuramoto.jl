##########################
#  Two species Kuramoto  #
##########################
using EcologicalNetworksDynamics
using Plots


function get_parameters_two_species_kuramoto(;
        ### Mortality rates for consumers
        dcons = .1,
        ### Handling time for consumers
        ht = 3.0,
        ### Attack rate
        ar = .7,
        # Preference for secondary resources
        B = .01,
        ## Growth rate of producers
        r = .3,
        ## Carrying capacity
        K = 1.0,
        ## Competition among producers
        α = .1
    )

    # Interaction matrix
    A = [
         0 0 0 0;
         0 0 0 0;
         1 1 0 0;
         1 1 0 0;
        ]

    # Preference matrix
    pref = [
            0 0 0 0;
            0 0 0 0;
            1 B 0 0;
            B 1 0 0;
           ]
    # Foodweb
    foodweb = FoodWeb(A)
    # Functional response
    classicresp = ClassicResponse(foodweb,
                                  # Hill exponent
                                  h = 1.0,
                                  # Consumer preference
                                  ω = pref,
                                  # Predator interference
                                  c = repeat([0], 4),
                                  # Handling time
                                  hₜ = ht,
                                  # Attack rates
                                  aᵣ = ar
                                 )
    # Biorates
    bio_rate = BioRates(foodweb,
                        d = [0, 0, dcons, dcons],
                        r = [r, r, 0, 0],
                        e = ones(4, 4),
                        x = 0.0,
                        y = 0.0
                       )
    # Environment
    env = Environment(foodweb, K = K)
    # Prod competition
    comp = [
            1 α 0 0;
            α 1 0 0;
            0 0 0 0;
            0 0 0 0;
           ]
    prod_comp = ProducerCompetition(foodweb, α = comp)

    # ModelParameters
    ModelParameters(foodweb,
                        functional_response = classicresp,
                        biorates = bio_rate,
                        environment = env,
                        producer_competition = prod_comp
                       );
end

# Trophic coupling
p = get_parameters_two_species_kuramoto(B = 0.01, α = 0)
m = simulate(p, [.1, .3, .1, .3],
             tmax = 2000,
             callback = nothing
            );
p1 = plot(m, idxs = [3, 4], legend = false, title = "Trophic coupling")
coefficient_of_variation(m, idxs = [3, 4])

# Resource coupling
p = get_parameters_two_species_kuramoto(B = 0.00, α = .1)
m = simulate(p, [.1, .3, .1, .3],
             tmax = 600,
             callback = nothing
            );
p2 = plot(m, idxs = [3, 4], legend = false, title = "Resource coupling")

# No Ressource  nor trophic coupling
p = get_parameters_two_species_kuramoto(B = 0.00, α = 0)
m = simulate(p, [.1, .3, .1, .3],
             tmax = 600,
             callback = nothing
            );
p3 = plot(m, idxs = [3, 4], legend = false, title = "No coupling")

myplot = plot(p1, p2, p3)

png(myplot, "kuramoto_two_species.png")
