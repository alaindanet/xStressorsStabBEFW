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

fw = FoodWeb(nichemodel,10, C = .3, Z = 100)
p = ModelParameters(fw, biorates = BioRates(fw; d = allometric_rate(fw, DefaultMortalityParams())))


#########
#  Sim  #
#########

df = DataFrame(Arrow.Table("scripts/param_comb_ct_S_h_d_allometric.arrow"))
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
            Z = 10, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 1000, last = 100, dt = 0.1, return_sol = true)

# test d
w = sim_int_mat(pm.A;
            ρ = 0.0, alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            K = 5.0,
            σₑ = 1,
            Z = 2, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 2000, last = 1000, dt = 0.1, return_sol = false, digits = 5)

using Distributed, ProgressMeter
timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (sim_id = p.sim_id, fw_id = p.fw_id,
                                productivity = p.K, h = p.h),
                               sim_int_mat(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = nothing,
                                           da = (ap = .4, ai = .4, ae = .4),
                                           σₑ = p.sigma, Z = p.Z,
                                           h = p.h, c = 0.0, K = p.K,
                                           dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
                                           max = 10000, last = 1000,
                                           K_alpha_corrected = true,
                                           dt = 0.1, gc_thre = .1,
                                           return_sol = false,
                                           digits = 5
                                          )
                              ),
                         param[repeat([1], 10)],
                         batch_size = 1
                        )
lala = DataFrame(sim);
lala[!, [:richness, :stab_com, :avg_cv_sp, :sync, :max_tlvl, :w_avg_tlvl, :avg_omnivory, :avg_int_strength]]


timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (sim_id = p.sim_id, fw_id = p.fw_id,
                                productivity = p.K, h = p.h),
                               begin
                                   fw = FoodWeb(p.A, Z = p.Z)
                                   S = richness(fw)
                                   funcrep = BioenergeticResponse(fw; h = p.h, c = 0)
                                   nprod = length(producers(fw))
                                   K = p.K * (1 + (.5 * (nprod - 1))) / nprod
                                   # Carrying capacity
                                   env = Environment(fw, K = K)
                                   d = allometric_rate(fw, AllometricParams(.4, .4, .4,-.25, -0.25, -.25))
                                   # generate the model parameters
                                   myp = ModelParameters(fw;
                                                         biorates = BioRates(fw; r = 1.0, d = d),
                                                         functional_response = funcrep,
                                                         environment = env,
                                                         producer_competition = ProducerCompetition(fw, αii = 1.0, αij = .5)
                                                        )

                                   B0 = rand(S)
                                   # Simulate
                                   m = simulate(myp, B0;
                                                alg = :adaptive,
                                                tmax = 10000,
                                                extinction_threshold = 1e-5,
                                                verbose = false
                                               );
                                   bm_cv = coefficient_of_variation(m; last = 1000)
                                   (
                                    richness = richness(m),
                                    max_tlvl = trophic_structure(m).max,
                                    stab_com = 1/bm_cv.community,
                                    avg_cv_sp = bm_cv.average_species,
                                    sync = bm_cv.synchrony
                                   )
                               end

                              ),
                         param[repeat([5], 10)],
                         batch_size = 1
                        )
DataFrame(sim)
fw = FoodWeb(nichemodel, 20, C = .22, check_cycle = true)
max_trophic_level(fw.A)
A= fw.A

w = sim_int_mat(A;
            ρ = 0.0, alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            K = 15.0,
            σₑ = 1.0,
            Z = 100, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 500, dt = .1, return_sol = false, digits = 5);
w
plot(w.species, tspan = (0, 500))
w = sim_int_mat(A;
            ρ = 0.0, alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            K = 20.0,
            σₑ = 1.5,
            Z = 1000, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 200, last = 100, dt = .1, return_sol = true, digits = 5);
EcologicalNetworksDynamics.check_last_extinction(w, idxs = 1:size(A,1), last = 100)
plot(w, idxs = 1:size(A, 1), tspan = (500, 1500))
plot(w, idxs = 1:size(A, 1))
w.t
tu = get_stab_fw(w; last = 1000)
tu.tlvl
connectance(A[tu.alive_species, tu.alive_species])

f = FoodWeb(nichemodel, 20, C = .22, check_cycle = true, Z = 100)

allometric_rate(f, AllometricParams(da.ap, da.ai, da.ae,-.25, -0.25, -.25))


# Check timesteps
w.tlvl
#

w = sim_int_mat(zeros(Int, 5, 5);
            ρ = 0.0, alpha_ij = 0.5,
            d = nothing,
            K = 5.0,
            σₑ = .5,
            Z = 100, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100, dt = 0.1, return_sol = false, digits = 5)

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

