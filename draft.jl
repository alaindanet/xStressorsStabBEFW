using Revise
using BEFWM2
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
using JSON3
using CSV
include("src/minmax.jl")
include("src/interaction_strength.jl")
include("src/stochastic_mortality_model.jl")
include("src/sim.jl")

ct = .1
S = 20

global fw
c = 0
try
    while c != ct
        fw = FoodWeb(nichemodel, 20, C = ct, Z = 4)
        c = round(connectance(fw), digits = 2)
    end
catch
   fw = missing
end
fw

#simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
ti = simCS(0.1, 20, 100, 2.0, 1.0, 1.0, 5; max = 5000, last = 1000, dt = 0.1, return_sol = false)
plot(ti, idxs = collect(1:1:length(get_parameters(ti).network.species)))
cv(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))
biomass(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))
BEFWM2.filter_sim(ti)
get_parameters(ti)




fw = FoodWeb([0 0 0; 1 0 0; 0 1 0])
p = ModelParameters(fw, functional_response = BioenergeticResponse(fw, h = 2, c = 1), env_stoch = EnvStoch(.3))
S = size(fw.A, 1)


stoch_starting_val = repeat([0], S)
u0 = [rand(S); stoch_starting_val]

# Make the stochastic matrix
corr_mat = zeros(S * 2, S * 2)
corr_mat[diagind(corr_mat)] .= 1.0
# Generate the Wiener Process
wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

prob = SDEProblem(
                  mydBdt!,
                  gen_stochastic_process,
                  u0,
                  [0, 50],
                  p,
                  noise = wiener
                 )
# Simulate
m = solve(prob;
      saveat = collect(0:1:50),
      dt = .1,
      adaptive = false
     )
cv(m, last = 10, idxs = collect(1:1:S))
biomass(m, last = 10, idxs = collect(1:1:S))
BEFWM2.filter_sim(m)
get_parameters(m)




ti = mysim(A, 1, 4, (h = 2, c = 1), 0, .3, max = 5000, last = 4000, dt = .1, corr_mat = vc, return_sol = true)

ti.stab_com
ti.tlvl

ti.int_strength
mean(ti.int_strength)
gini(([i for i in ti.int_strength if i != 0]))
ti.max_int
ty = [i for i in ti.max_int if i != 0]
minimum(ty)
gini(vec(ty))
gini(ty)
sum(ti.tlvl .* ti.bm_sp/sum(ti.bm_sp))
sum([1, 3] .* [3, 1]/sum([3,1]))

df = DataFrame([ti, ti])
typeof(df)

CSV.write("csv_test", df)
write(write("df_test", String), df)
tu = CSV.read("csv_test", DataFrame)
tu[1, :bm_sp]


DataFrame(tu.columns)

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
