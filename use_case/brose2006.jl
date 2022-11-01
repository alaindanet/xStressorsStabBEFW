
reps = 5
S = [20, 30, 40]
# list to store networks
global networks = []
# monitoring variable 
for s in 1:length(S)
    l = 0 
    # while loop
    while l < reps
        # generate network
        fb = FoodWeb(nichemodel, S[s], C = 0.15)
        # convert the UnipartiteNetwork object into a matrix of 1s and 0s
        A = Int.(fb.A)
        # calculate connectance
        co = sum(A)/(size(A,1)^2)
        # ensure that connectance = 0.15
        if co == 0.15
            push!(networks, A)
            l = l + 1
            # save network is co = 0.15
        end
    end
end

# three structural food web models, three functional responses, two metabolic
# categories, three levels of diversity and eight body mass ratios whose logs are
# evenly spaced from 10)2 to 105
#
# Functional response
fr_type = [(h = 1.0, c = 0.0, name = "Type II"),
           (h = 2.0, c = 0.0, name = "Type III"),
           (h = 1.0, c = 1.0, name = "Type II with predator interference")]
Z = range(10^-2,10^5, length = 5)



# Storage
df = []
#z, f, h = Z[1], fr_type[1], 1
for h in 1:reps
    A = networks[h]
    for f in fr_type 
        for z in Z 

            fb = FoodWeb(A, Z = z)
            # create model parameters
            
            # Prod
            prod_mask = [sum(fb.A[i,:]) == 0 for i in 1:size(fb.A, 1)]
            
            ## BioenergeticResponse
            # Equal comsumption rate among consumers for a given prey:
            ω_t = mapslices(x -> x / sum(x),fb.A, dims = 2)

            ## Biorates
            ### Conversion efficiency of predators on prey
            nsp = size(fb.A, 1)
            E = spzeros(nsp, nsp)
            for i in 1:nsp
                for j in 1:nsp
                    if fb.A[i,j] == 0# Si no interaction, 0
                        E[i, j] = 0
                    elseif sum(fb.A[j,:]) == 0 #if the prey is a plant
                        E[i,j] = .45 # Predator conversion from plant
                    else
                        E[i,j] = .85 # Predator conversion from non-plant
                    end
                end
            end
            ### Growth rate = 1.0 only for producers 
            r = [if prod_mask[i] == true 1.0 else 0.0 end for i in 1:size(prod_mask, 1)] 
            #
            p = ModelParameters(fb,
                            functional_response =
                            BioenergeticResponse(fb,
                                                 h = f.h,
                                                 c = f.c,
                                                 # Equal comsumption rates among
                                                 # consumers for a given prey:
                                                 ω = ω_t
                                                                      ),
                            environment = Environment(fb, K = 1.0),
                            biorates = BioRates(fb; e = E, r = r)
                           );
        # assign biomasses between 0.05 and 1
        bm = rand(5:100, size(A,1)) / 100
        # simulate
        out = simulate(p, bm, tmax=2000,
                       callback = CallbackSet(
                                              BEFWM2.PositiveDomain(),
                                              # TerminateSteadyState(1e-6, 1e-4),
                                              BEFWM2.ExtinguishSpecies(10^-30,false),
                                             )
                      )

        # dummy naming variables
        # save `out` as a JLD2 object using the @save macro:
        # @save "out_objects/model_output, network = $h, alpha = $α_num, K = $K_num.jld2" out

        # calculate output metrics
        diversity = species_persistence(out, threshold = 10^-30, last = 1)
        stability = foodweb_cv(out, threshold = 10^-30, last = 1000)
        biomass = total_biomass(out, last = 1000)

        # push to df
        rich = size(A, 1)
        push!(df, (fr = f.name, Z = z, S = rich, network = "$h", persistence = diversity, cv = stability, bm = biomass))

        # print some stuff - see how the simulation is progressing
        println(("f = $f", "Z = $z", "S = $rich", "network = $h"))
        end
    end
end

ti = DataFrame(df)

plot(ti[:, "persistence"], ti[:, "cv"])
