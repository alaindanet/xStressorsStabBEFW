using Revise
using BEFWM2
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
include("src/minmax.jl")


# Try the solver option with 
function connectance(mat)
    sum(mat) / size(mat,1)^2
end

foodweb = FoodWeb(nichemodel, 10,C = .14)
while connectance(foodweb.A) != .14 
    global foodweb = FoodWeb(nichemodel, 10,C = .14)
end
# Generally, reltol is the relative accuracy while abstol is the accuracy when u
# is near zero. These tolerances are local tolerances and thus are not global
# guarantees. However, a good rule of thumb is that the total solution
# accuracy is 1-2 digits less than the relative tolerances. Thus for the
# defaults abstol=1e-6 and reltol=1e-3, you can expect a global accuracy
# of about 1-2 digits.
q0 = 0.0:.025:.5
params = ModelParameters(foodweb, functional_response = BioenergeticResponse(foodweb, h = 0.025 + 1))
sol = simulate(
               params,
               rand(size(foodweb.A, 1));# repeat([.1], size(foodweb.A, 1)),
               extinction_threshold = 1e-10, 
               callback = CallbackSet(
                                      BEFWM2.PositiveDomain(),
                                      TerminateSteadyState(1e-6, 1e-4),
                                      BEFWM2.ExtinguishSpecies(1e-10, false),
                                     ),
               reltol=1e-8, abstol=1e-8
              )
using DifferentialEquations
using OrdinaryDiffEq


function connectance(mat)
    sum(mat) / size(mat,1)^2
end

Random.seed!(1234)
foodweb = FoodWeb(nichemodel, 10,C = .14)
# I would like a fw with c = 0.14 plz
while connectance(foodweb.A) != .14 
    global foodweb = FoodWeb(nichemodel, 10,C = .14)
end
params = ModelParameters(foodweb, functional_response = BioenergeticResponse(foodweb, h = 0.025 + 1))
B = rand(size(foodweb.A, 1))

# Solvers
non_stiff_solver = [Vern9(), OwrenZen3(), BS5(), Tsit5()]
stiff_solver = [Rosenbrock23(), Rodas4P(), TRBDF2()] 

for a in [non_stiff_solver; stiff_solver]
    println("solver is: $(a)")
    try
        sol = simulate(
                       params,
                       B,
                       alg = a,
                       tmax = 3000,
                       callback = CallbackSet(
                                              BEFWM2.PositiveDomain(),
                                              #TerminateSteadyState(1e-6, 1e-4),
                                              BEFWM2.ExtinguishSpecies(1e-10, false),
                                             ),
                      )
    catch err
        println("Error with: $(err) \n")
    end
end

# Regarding the domain error with Tsit5(): DomainError(-4.2275272782171454e-7, "Exponentiation yielding a complex result...
# A negative number connot be powered by a decimal
a = -1
(a)^1.025
(-1)^1.025
# But work with precedence rules (power before negative)
-1^1.025

Iₖ = ωₖ * x * y / B₀

