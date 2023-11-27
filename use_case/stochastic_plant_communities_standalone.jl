using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Distributions
using Statistics
using Plots


# dBdt of plants and Stochastic process 
function f(du, u, p, t)
    # Logistic growth and natural mortality
    for i in 1:S
        du[i] = r * u[i] * (1 - u[i] / K) - d * u[i] * exp(u[S+i])
    end
    # Stochastic variation in mortality as a *** process (Wiener process that stays around
    # the mean)
    for i in S+1:2*S
        du[i] = θ * (0 - u[i])
    end
end

# Stochastic brownian process of standard deviation sigma
function g(du, u, p, t)
    du .= 0.0
    du[S+1:2*S] .= σ
end

# Simulation parameters 
S = 2
r = 1.0
K = 10.0 / S
d = .4
θ = 1
σ = .6
ρ = 0.0

# Build the correlated wiener process
stoch_starting_val = repeat([0], S)
# Make the stochastic matrix
corr_mat = zeros(S * 2, S * 2)
corr_mat .= ρ
corr_mat[diagind(corr_mat)] .= 1.0
tmax = 50000
tspan = (0.0, tmax)
heston_noise = CorrelatedWienerProcess!(corr_mat,
                                        tspan[1],
                                        zeros(size(corr_mat, 1)),
                                        zeros(size(corr_mat, 1))
                                       )

# Starting values and problem definition
u0 = [rand(S); stoch_starting_val]
pb = SDEProblem(f, g, u0, tspan, noise = heston_noise)
# Solve the dynamic 
m = solve(pb,
      saveat = collect(0:1:tmax),
      dt = .1,
      adaptive = false)
# See the stochastic process
plot(1:tmax, transpose(m[S+1:2*S,1:tmax]))
# See the biomass
plot(1:tmax, transpose(m[1:S,1:tmax]))
