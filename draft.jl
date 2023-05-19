using Revise
using EcologicalNetworksDynamics
using EcologicalNetworksDynamics: richness
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Distributions
using Statistics
using DataFrames
using Plots
using Debugger
using CSV
using Arrow
include("src/minmax.jl")
include("src/interaction_strength.jl")
#include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/plot.jl")
include("src/get_modules.jl")

files = readdir("res/simCSZenrich2")
df = DataFrame(Arrow.Table(join(["res/simCSZenrich2/", files[1]])))

select(df, [:species, :stoch])
for i in files[isfile.([join(["res/simCSZenrich2/", i]) for i in files])]
    df = DataFrame(Arrow.Table(join(["res/simCSZenrich2/", i])))
    select!(df, Not([:cv_sp]))
    Arrow.write(join(["res/simCSZenrich2/simCSZenrich2/", i]), df, ntasks = 2)
end


arrow = Arrow.Table("res/simCSZenrich2/simCSZ_4001_8000.arrow")
df = DataFrame(arrow)
filter(:rho => x -> x == 1.0, df)[3, [:cv_sp, :bm_sp]]
filter(:cv_sp => x -> ismissing(x), df)
filter(:cv_sp => x -> !ismissing(x) , df)
df.cv_sp

toy = select(df, [:rep, :productivity, :richness, :ct, :Z, :rho, :env_stoch, :species, :stoch])
Arrow.write("res/simCSZenrich2/simCSZenrich2/sim_toy_compensatory_dyn.arrow", toy)

#########
#  Sim  #
#########

#simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)
using Distributions

ti = simCS(.4, 5;
           d = .2,
           Z = 1, h = 2.0, ρ = -1,
           σₑ = 1.0, K = .5, max = 500,
           last = 100, dt = 0.1, return_sol = false
          )
map(typeof, ti)
convert.(Float64, ti.omnivory)


nb_alive_species = try
    length(trophic_structure(ti, last = 1000).alive_species)
catch
    missing
end

plot(ti, idxs = collect(1:1:length(get_parameters(ti).network.species)))
cv(ti, last = 100, idxs = collect(1:1:length(get_parameters(ti).network.species)))
biomass(ti, last = 10, idxs = collect(1:1:length(get_parameters(ti).network.species)))



# fw = FoodWeb([0 0 0; 1 0 0; 0 1 0], Z = 100)
fw = FoodWeb([0 0 0; 0 0 0; 0 0 0], Z = 100)
fw = FoodWeb([0 0; 0 0], Z = 100)
p = ModelParameters(fw,
                    functional_response = BioenergeticResponse(fw, h = 2, c = 1),
                    producer_competition = ProducerCompetition(fw; αij = .5),
                    env_stoch = EnvStoch(.5),
                    biorates = BioRates(fw; d = 0.05)
                   )
S = size(fw.A, 1)

stoch_starting_val = repeat([0], S)
u0 = [0.0, 1]
m = simulate(p, u0;
         rho = 1,
         dt = .5,
         tmax = 50,
         extinction_threshold = 1e-5,
         verbose = false
        );

get_stab_fw(m; last = 10)
cv(m, last = 100, idxs = collect(1:1:S))
EcologicalNetworksDynamics.synchrony(transpose(m[S+1:1:2*S, end-(100-1):end]))

biomass(m, last = 100, idxs = collect(1:1:S))
plot(m, idxs = collect(1:1:S))
plot(m, idxs = collect(S+1:1:2 * S))

sim_int_mat([0 0 0; 0 0 0; 1 1 0];
            ρ = 1.0, alpha_ij = 0,
            d = 0.1,
            σₑ = .5, Z = 100, h = 2.0, c = 1.0, K = 1.0,
            fun = stoch_d_dBdt!,
            max = 500, last = 100, dt = 0.1, return_sol = false)



###########
#  Motif  #
###########

fw = FoodWeb(nichemodel, 20, C = .05)
fw = FoodWeb(nichemodel, 80, C = .05)
EcologicalNetworks.find_motif(fw.A, unipartitemotifs().S1)
mot = find_motif(UnipartiteNetwork(fw.A), unipartitemotifs().S1) |> length
map(x -> find_motif(UnipartiteNetwork(fw.A), x) |> length, unipartitemotifs())

# Diamond
map(x -> find_motif(UnipartiteNetwork([0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0] .> 0), x) |> length, unipartitemotifs())

##########
#  Plot  #
##########

webplot(get_fw_modules()[10]; consasrow = true)
webplot([0 0 0; 1 1 0; 1 1 0]; consasrow = true)
ti = map(x -> webplot(x; consasrow = true), get_fw_modules())
length(ti)
plot(ti..., layout = (4, 4))
plot(ti[1], ti[2], ti[3], layout = (1, 3))
plot((ti[i] for i in 1:length(ti))..., layout = (4, 3))

