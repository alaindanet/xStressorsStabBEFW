using Revise
using EcologicalNetworksDynamics
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

#########
#  Sim  #
#########

df = DataFrame(Arrow.Table("scripts/param_comb_ct_S_h.arrow"))
reshape_array(vec) = reshape(vec, (
                                   round(Int, sqrt(length(vec))),
                                   round(Int, sqrt(length(vec)))
                                  )
                            )
df[!, :A] = map(x -> reshape_array(x), df[:, :A])

param = NamedTuple.(eachrow(df))

pm = sample(param)
println("Running warmup: K = $(pm.K), σₑ = $(pm.sigma), ρ = $(pm.rho)")

warmup = sim_int_mat(pm.A;
            ρ = 1.0, alpha_ij = 0,
            d = 0.0,
            σₑ = .1, Z = 1000, h = 2.0, c = 0.0, K = 1.0,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 50, last = 10, dt = 0.1, return_sol = false, digits = 5)
println("$(warmup)")

w = sim_int_mat(pm.A;
            ρ = 0.0, alpha_ij = 0,
            d = 0.0,
            K = 1.0,
            σₑ = .1,
            Z = 1000, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 1000, last = 100, dt = 0.1, return_sol = true)
tu = get_stab_fw(w, digits = 5)
tu.alive_species
varinfo(r"tu")
[3,2,1][findall([false,true,false])[1]]

[3,2,1][]
[3,2,1][[false,true,false]]

tui = get_stab_fw(w)
varinfo(r"tui")

w = sim_int_mat(pm.A;
            ρ = 1.0, alpha_ij = 0,
            d = 0.0,
            K = 5.0,
            σₑ = .5,
            Z = 100, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 1000, last = 100, dt = 0.1, return_sol = false, digits = 5)
# Check timesteps

w.tlvl
w.alive_species
w.bm_sp
w.cv_sp
w.met_loss

using ProgressMeter, Distributed
sim = @showprogress pmap(p ->
                         merge(
                               (fw_id = p.fw_id, productivity = p.K, h = p.h),
                               sim_int_mat(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = 0.0,
                                           σₑ = p.sigma, Z = p.Z,
                                           h = p.h, c = 0.0, K = p.K,
                                           dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
                                           max = 5000, last = 1000,
                                           K_alpha_corrected = true,
                                           dt = 0.1, gc_thre = .2,
                                           return_sol = false
                                          )
                              ),
                         param[collect(1:5:100)],
                         batch_size = 2
                        )
df_res = DataFrame(sim)

failed_sim = @showprogress pmap(p ->
                         merge(
                               (fw_id = p.fw_id, productivity = p.K, h = p.h),
                               sim_int_mat(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = 0.0,
                                           σₑ = p.sigma, Z = p.Z,
                                           h = p.h, c = 0.0, K = p.K,
                                           dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
                                           max = 5000, last = 1000,
                                           K_alpha_corrected = true,
                                           dt = 0.5, gc_thre = .2,
                                           return_sol = false
                                          )
                              ),
                         param[collect(1:4:10)],
                         batch_size = 2
                        )

df_failed = DataFrame(failed_sim)


#############
#  testing  #
#############

sim_int_mat([0 0 0; 0 0 0; 1 1 0];
            ρ = 1.0, alpha_ij = 0,
            d = 0.0,
            σₑ = .5, Z = 1000, h = 2.0, c = 0.0, K = 1.0,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 1000, last = 500, dt = 0.1, return_sol = false)

###########
#  Motif  #
###########

fw = FoodWeb(nichemodel, 20, C = .2, Z = 10)
p = ModelParameters(fw)
p.biorates.x
trophic_levels(fw)
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

