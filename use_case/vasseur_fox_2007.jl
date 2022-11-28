using Revise
using BEFWM2
using Plots
using LinearAlgebra
using DifferentialEquations
using DataFrames
using CSV
include("../src/stochastic_mortality_model.jl")

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

bioener = BioenergeticResponse(foodweb,
                               h = 1,
                               # Half saturation-constant
                               B0 = [0, 0.16129, .9, .5],
                               # Consumer preference
                               ω = myω,
                               # Predator interference
                               c = myc
                              )

# Metabolic rates and maximum ingestion rates (merged in J parameter in Vasseur
# & Fox)
J = [0, 0.8036, 0.7, .4]
## From McCann (1998):
x = [0, .40, .20, .08]
y = [0, 2.009, 3.50, 5.0]
x .* y .- J

biorate = BioRates(foodweb,
        r = [1.0, 0, 0, 0],
        e = mye,
        x = x,
        y = y
       )




# define the starting values
S = BEFWM2.richness(foodweb)
stoch_starting_val = [0; 0; 0; 0]
u0 = [rand(S); stoch_starting_val]
tspan = (0, 50000)

# Simulation
#
# Along correlation between environmental stochasticity and strength of
# stochasticity
ρ = collect(-1.0:0.1:1.0)
σₑ = collect(0.0:.05:.6)
rep = collect(1:5)
tspan = [0, 50000]
tcol = collect(0:1:tspan[2])
dt = 0.1
nb_ts_for_cv = 25000

# Build the covariance matrix of the size of two times the number of species
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc

# define the starting values
stoch_starting_val = [0; 0; 0; 0]
corr_mat = vc

out = []
for p in ρ
    for s in σₑ
        for r in rep
            println("Start ρ = $p, σₑ = $s, rep = $r")

            B0 = [rand(S); stoch_starting_val]
            params_tmp = ModelParameters(foodweb,
                                         functional_response = bioener,
                                         biorates = biorate,
                                         environment = Environment(foodweb, K = 1.0),
                                         env_stoch = EnvStoch(s)
                                        )
            # Correlation between consumer
            corr_mat[S+2, S+3] = p
            corr_mat[S+3, S+2] = p

            # Generate the Wiener Process
            wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

            prob = SDEProblem(
                              mydBdt!,
                              stochastic_process,
                              u0,
                              tspan,
                              params_tmp,
                              noise = wiener
                             )
            sol = solve(prob;
                        saveat = tcol,
                        dt = dt,
                        adaptive = false
                       )
            # TODO: Add extinction check
            if length(sol.t) == length(tcol)

                cv = foodweb_cv(sol, last = nb_ts_for_cv, idxs = [1, 2, 3, 4])
                sync_cons = foodweb_cv(sol, last = nb_ts_for_cv, idxs = [2, 3]).synchrony

                push!(out, (
                            ρ = p, σₑ = s, rep = r, #sol = sol,
                            stab_com = 1 / cv.cv_com, avg_cv_sp = cv.avg_cv_sp, sync = cv.synchrony,
                            stab_cons1 = 1 / cv.cv_sp[2], stab_cons2 = 1 / cv.cv_sp[3], stab_pred = 1 / cv.cv_sp[4],
                            sync_cons = sync_cons
                           )
                     )
            else
                push!(out, (
                            ρ = p, σₑ = s, rep = r, #sol = sol,
                            stab_com = nothing, avg_cv_sp = nothing, sync = nothing,
                            stab_cons1 = nothing, stab_cons2 = nothing, stab_pred = nothing,
                            sync_cons = nothing
                           )
                     )

            end
        end
    end
end

df = DataFrame(out)

CSV.write("vasseur_fox_res.csv", df)
