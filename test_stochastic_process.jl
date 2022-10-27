
using DifferentialEquations
using DiffEqFinancial
HestonProblem(μ,κ,Θ,σ,ρ,u0,tspan)
ti = HestonProblem(2, 1, 1, 1, .5, 0, 1)

μ = 1.0
σ = 2.0
W = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)
# ...
# Define f,g,u0,tspan for a SDEProblem
# ...
prob = SDEProblem(f,g,u0,tspan,noise=W)

Θ = 1.0
μ = 0.1
σ = 2.0
W = OrnsteinUhlenbeckProcess(Θ, μ, σ, 0.0, 1.0, 1.0)
prob = NoiseProblem(W,(0.0,100))
sol = solve(prob;dt=0.1)
plot(sol.u)

# Wiener process WienerProcess(t0,W0,Z0=nothing;kwargs...)
W = WienerProcess(0.0, 0.0)
prob = NoiseProblem(W,(0.0,1000))
sol = solve(prob; dt=1.0)
plot(sol.u)

using LinearAlgebra, Random
b = LinearAlgebra.svd([1 2; 3 4])
b.S, b.U

γ = LinearAlgebra.svd([1 -1; -1 1])
γ.S, γ.U
A = γ.U * Diagonal(sqrt.(γ.S))


gi = A * sqrt.(abs(4.3)) * randn(Random.MersenneTwister(1234), Float64, (2, 100))
plot(transpose(gi))
