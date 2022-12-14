using Revise
using BEFWM2
using Plots
using LinearAlgebra
using DifferentialEquations
using DataFrames
using CSV
include("src/stochastic_mortality_model.jl")

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

# Build the covariance matrix of the size of two times the number of species
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc

# define the starting values

corr_mat = vc

function simulate_stoch(foodweb, bioener, biorate, s, p, corr_mat; tmax = 1000, dt = .1)

    S = BEFWM2.richness(foodweb)
    stoch_starting_val = [0; 0; 0; 0]
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
                      B0,
                      [0, tmax],
                      params_tmp,
                      noise = wiener
                     )
    solve(prob;
          saveat = collect(0:1:tmax),
          dt = dt,
          adaptive = false
         )
end

function get_ts(solution; last = 400, rho = missing, sigma = missing)
    ts = transpose(solution[1:4, end-(last - 1):end])
    df = DataFrame(ts, ["resource", "cons1", "cons2", "pred"])
    df[!, :time] = 1:last
    df[!, :sigma] .= sigma
    df[!, :rho] .= rho
    df
end

sol.prob.p.env_stoch.σₑ

sol = simulate_stoch(foodweb, bioener, biorate, 0, 0, vc, tmax = 500)
df_no_stoch = get_ts(sol, last = 400, rho = 0, sigma = 0)
plot(sol, idxs = [1, 2, 3, 4])
foodweb_cv(sol, last=400, idxs = [1, 2, 3, 4])

sol_stoch_cor0 = simulate_stoch(foodweb, bioener, biorate, .3, 0, vc, tmax = 500)
df_stoch_cor0 = get_ts(sol_stoch_cor0, last = 400, rho = 0, sigma = 0.3)
plot(sol, idxs = [1, 2, 3, 4])
foodweb_cv(sol, last=400, idxs = [1, 2, 3, 4])

sol_stoch_cor1 = simulate_stoch(foodweb, bioener, biorate, .3, 1, vc, tmax = 500)
df_stoch_cor1 = get_ts(sol_stoch_cor1, last = 400, rho = 1, sigma = 0.3)
plot(sol, idxs = [1, 2, 3, 4])
foodweb_cv(sol, last=400, idxs = [1, 2, 3, 4])


df = [df_no_stoch; df_stoch_cor0; df_stoch_cor1]
CSV.write("./use_case/vasseur_fox_sim_ts_ex.csv", df)
