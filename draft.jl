using Revise
using BEFWM2, SparseArrays
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
include("src/minmax.jl")


fw = FoodWeb(nichemodel, 20, C = .15, Z = 1000)
# fw = FoodWeb(nichemodel, 5, C = .1)
p = ModelParameters(fw)

#Iₖ = ωₖ * x * y / B₀

ti = simulate(p, rand(size(p.network.species, 1)))

foodweb = FoodWeb([0 0; 1 0]); # create the foodweb
params = ModelParameters(foodweb); # generate the parameters
B = [0.5, 0.5]; # set initial biomass
solution = simulate(params, B); # run simulation
typeof(solution)
params.functional_response.B0
params.biorates.x

emp_int_str = empirical_interaction_strength(bm, p)
gini(vec(emp_int_str))
