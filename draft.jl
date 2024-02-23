using EcologicalNetworksDynamics: simulate_deter
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
using Distributed, ProgressMeter
include("src/minmax.jl")
include("src/interaction_strength.jl")
#include("src/stochastic_mortality_model.jl")
include("src/sim.jl")
include("src/plot.jl")
include("src/get_modules.jl")

fw = FoodWeb(nichemodel, 10, C = .1, Z = 1)
bioener = BioenergeticResponse(fw,
                               h = 3,
                               # Predator interference
                               c = 0
                              )
p = ModelParameters(fw,
                    functional_response = bioener,
                    env_stoch = EnvStoch(.3))
tmax = 500
m = simulate(p, rand(richness(fw));
         rho = 0,
         dt = .1,
         tmax = tmax
        )

species_persistence(m)

pm = get_parameters(m)
alive_species = trophic_structure(m, last = 1).alive_species
check_disconnected_species(pm.network.A, alive_species)

ti = [
 0 0 0 0;
 0 0 0 0;
 1 0 0 0;
 0 1 0 0
]
pu = check_disconnected_species(ti, [1, 2, 3, 4])


bm = [1, 1, 3, 4]
alive = [1, 3, 4]
kill_disconnected_species(ti, alive_species = [1, 2, 3, 4], bm = bm) == [1, 1, 3, 4]
kill_disconnected_species(ti, alive_species = [], bm = bm) == [0, 0, 0, 0]
kill_disconnected_species(ti, alive_species = [1, 3, 4], bm = bm) == [1, 0, 3, 0]
kill_disconnected_species(ti, alive_species = [1, 2, 4], bm = bm) == [0, 1, 0, 4]
kill_disconnected_species(ti, alive_species = [1, 2, 4], bm = bm) == [0, 1, 0, 4]

to = [1, 2, 3, 4]
tu = to[.![false, false, false, false]]
bm_to_set_to_zero = 1:length(bm) .∉ [tu]
bm[bm_to_set_to_zero] .= 0
1:length(bm)

idxs = alive_species
ti = filter_model_parameters(p, idxs = alive_species)
ti.environment.K

coefficient_of_variation(m, last = 100)

simulate(p, rand(richness(fw)), saveat = (0:1:tmax))

m = (t = 0, b = 0)
m[:t]
dictionary = Dict(1 => 77, 2 => 66, 3 => 1)
keys(dictionary)
dictionary[:1]

param = DataFrame(Arrow.Table("scripts/param_comb_ct_S_h_d3.arrow"))

filter([:S, :sigma, :h] => (rich, s, h) -> rich == 40 && s == .6 && h == 2, param)

#toy_param = param[[2, 30, 58, 362, 390, 742, 750, 766, 1082, 1098],:]
toy_param = param[[1, 3, 4, 6, 7, 10, 11, 13, 14],:]

# Reshape interaction matrix
reshape_array(vec) = reshape(vec, (
                                   round(Int, sqrt(length(vec))),
                                   round(Int, sqrt(length(vec)))
                                  )
                            )
toy_param[!, :A] = map(x -> reshape_array(x), toy_param[!, :A])

# Make a tuple vector
toy_param = NamedTuple.(eachrow(toy_param))

include("src/sim.jl")
timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (sim_id = p.sim_id, fw_id = p.fw_id, h = p.h),
                               sim_int_mat_check_disconnected(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = nothing,
                                           da = (ap = .4, ai = .4, ae = .4),
                                           σₑ = p.sigma, Z = p.Z,
                                           h = p.h, c = 0.0, K = 10,
                                           dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                                           max = 5000, last = 500,
                                           K_alpha_corrected = true,
                                           dt = 0.1, gc_thre = .1,
                                           dt_rescue = .05,
                                           return_sol = false,
                                           re_run = false,
                                           digits = 5
                                          )
                              ),
                         toy_param,
                         batch_size = 100
                        )
sim_df = DataFrame(sim)
select(sim_df, [:max_tlvl, :richness])
names(sim_df)


ti = map(x -> remove_disconnected_species(x.omega, x.alive_species),sim)
tu = FoodWeb(ti[4])
tu = FoodWeb(ti[9])
tu.A
EcologicalNetworksDynamics.is_connected(tu)
tu.A
tu.metabolic_class


ti = sim_int_mat_check_disconnected(toy_param[6].A;
            ρ = 0,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 3.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 500, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = false,
            re_run = true,
            digits = 5
           )

ti = sim_int_mat_check_disconnected(toy_param[1].A;
            ρ = 0,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 3.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 500, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = false,
            re_run = true,
            digits = 5
           )
ti.species

p = get_parameters(m)
p.network.A
plot(m.species)
FoodWeb(m.omega .> 0)
FoodWeb(m.int_strength .> 0)
m.int_strength
ω

p = ModelParameters(tu)
m = simulate_deter(p, rand(richness(tu)); tmax = 500)
m.u
ty = m[:, end]
ty[ty .> 0]

plot(m)

tmp_A = sim[4].omega .> 0
alive_species = sim[4].alive_species
cons = sum.(eachrow(tmp_A)) .!= 0
prod = sum.(eachrow(tmp_A)) .== 0
living_A = tmp_A[alive_species, alive_species]
alive_cons = cons[alive_species]
alive_prod = prod[alive_species]

species_with_no_pred = sum.(eachcol(living_A)) .== 0
disconnected_prod = alive_prod .&& species_with_no_pred
disconnected_cons = sum.(eachrow(living_A)) .== 0 .&& alive_cons
to_keep = .!(disconnected_prod .|| disconnected_cons)

living_A[to_keep, to_keep]
tmp_cons = sum.(eachrow(tmp_A)) .!= 0
sp_with_no_pred = sum.(eachcol(living_A)) .== 0
disconnected_prod = alive_prod .&& sp_with_no_pred
disconnected_cons = sum.(eachrow(living_A)) .== 0 .&& alive_cons

FoodWeb(tmp_A)
FoodWeb(tmp_A[.!( prod .&& sp_with_no_pred), .!( prod .&& sp_with_no_pred)])
select(sim_df, [:richness, :w_avg_tlvl, :max_tlvl])
test_param = toy_param[6]
ti = sim_int_mat(test_param.A;
            ρ = test_param.rho,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = test_param.sigma, Z = test_param.Z,
            h = test_param.h, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 5000, last = 500,
            K_alpha_corrected = true,
            extinction_threshold = 1e-5,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = true,
            re_run = true,
            digits = 5
           )

ti.retcode == DifferentialEquations.ReturnCode.Success
te = get_time_series(ti; last = 5000, idxs = nothing, digits = 10)
plot(1:5000, std.(eachrow(te.stoch)))
plot(0:5000, minimum.(eachrow(te.species)))
minimum.(eachcol(te.species))

tmax = 20000
S = 2
A = zeros(Int64, S, S)
ti = sim_int_mat(A;
            ρ = 0.0,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .6, Z = 10,
            h = 2.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = tmax, last = 500,
            K_alpha_corrected = true,
            dt = 0.5, gc_thre = .1,
            return_sol = true,
            re_run = false,
            dt_rescue = nothing,
            digits = 5
           )
ti.retcode == DifferentialEquations.ReturnCode.Unstable
ti = sim_int_mat(A;
            ρ = 0.0,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .6, Z = 10,
            h = 2.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 500,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            return_sol = true,
            re_run = false,
            digits = 5
           )
ti.retcode === [:Success]
ti.retcode == 1
ti.retcode.Success
names(ti.retcode)
te = get_time_series(ti; last = tmax, idxs = nothing, digits = 10)
plot(1:tmax, te.stoch)
plot(1:tmax, transpose(ti[S+1:2*S,1:tmax]))


tmax = 5000
S = 10; Z = 10
A = FoodWeb(nichemodel, 10; C = .3, Z = Z, quiet = true).A
ti = sim_int_mat(A;
            ρ = 0.0,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .6, Z = 10,
            h = 2.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = tmax, last = 500,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            return_sol = false,
            re_run = true,
            digits = 5
           )
plot(ti.species)
plot(ti.stoch)

A = toy_param[1].A
Z = toy_param[1].Z
alpha_ij = .5
K = 10
K_alpha_corrected = true
da = (ap = .4, ai = .4, ae = .4)
d = nothing
σₑ = .6
ρ = 0
dt = .1
r = 1
tmax = 500

fw = FoodWeb(A, Z = Z, quiet = true)
S = richness(fw)
C = connectance(fw)

# Parameters of the functional response
funcrep = BioenergeticResponse(fw; h = 2.0, c = 0.0)
# Carrying capacity
env = Environment(fw,
                  K = if K_alpha_corrected
                      nprod = length(producers(fw))
                      # Loreau & de Mazancourt (2008), Ives et al. (1999)
                      K * (1 + (alpha_ij * (nprod - 1))) / nprod
                  else
                      K
                  end
                 )
if isnothing(d)
    if !isnothing(da)
        d = allometric_rate(fw, AllometricParams(da.ap, da.ai, da.ae,-.25, -0.25, -.25))
    else
        d = allometric_rate(fw, DefaultMortalityParams())
    end
end
# generate the model parameters
p = ModelParameters(fw;
                    biorates = BioRates(fw; r = r, d = d),
                    functional_response = funcrep,
                    environment = env,
                    producer_competition = ProducerCompetition(fw, αii = 1.0, αij = alpha_ij),
                    env_stoch = EnvStoch(σₑ)
                   )

ω = p.functional_response.ω
B0 = rand(S)
# Simulate
m = simulate(p, B0;
             rho = ρ,
             dt = dt,
             tmax = tmax,
             extinction_threshold = .1,
             diff_code_data = (EcologicalNetworksDynamics.stoch_d_dBdt!, p),
             verbose = true
            );
plot(1:tmax, transpose(m[S+1:2*S,1:tmax]))
plot(1:tmax, transpose(m[1:S,1:tmax]))


# Define ODE problem and solve
t0 = 0.0
diff_code_data = (EcologicalNetworksDynamics.stoch_d_dBdt!, p)

timespan = (t0, tmax)
timesteps = collect(t0:dt:tmax)
code, data = diff_code_data

# Work around julia's world count:
# `generate_dbdt` only produces anonymous code,
# so the generated functions cannot be overriden.
# As such, and in principle, the 'latest' function is unambiguous.
if !isa(code, Function)
    message = "The given specialized code is not a `Function` but `$(typeof(code))`."
    if isa(code, Expr)
        message *= " Did you forget to `eval()`uate it before passing it to `simulate()`?"
    end
    throw(ArgumentError(message))
end
fun = (args...) -> Base.invokelatest(code, args...)
extinct_sp = Dict(i => 0.0 for (i, b) in enumerate(B0) if b == 0.0)
myp = (params = data, extinct_sp = extinct_sp)

stoch_starting_val = repeat([0], S)
u0 = [rand(S); stoch_starting_val]
# Make the stochastic matrix
corr_mat = zeros(S * 2, S * 2)
corr_mat .= rho
corr_mat[diagind(corr_mat)] .= 1.0
wiener = CorrelatedWienerProcess(corr_mat, t0, zeros(size(corr_mat, 1)))

problem = SDEProblem(
                     code,#stoch_d_dBdt!,
                     EcologicalNetworksDynamics.gen_stochastic_process,
                     u0,
                     timespan,
                     myp,
                     noise = wiener
                    )
# Try to simulate
m = solve(problem;
      saveat = collect(t0:1:tmax),
      dt = dt,
      adaptive = false,
      # callback = callback,
      isoutofdomain = (u, p, t) -> any(x -> x < 0, u)
     )
# Works 
plot(1:tmax, transpose(m[S+1:2*S,1:tmax]))
plot(1:tmax, transpose(m[1:S,1:tmax]))

# Does not work:
m = solve(problem;
      saveat = collect(t0:1:tmax),
      dt = dt,
      adaptive = false,
      callback = CallbackSet(
        #TerminateSteadyState(1e-6, 1e-4),
        ExtinctionCallback(1e-5, p, true),
       ),
      isoutofdomain = (u, p, t) -> any(x -> x < 0, u[1:S])
     )
plot(1:tmax, transpose(m[S+1:2*S,1:tmax]))
plot(1:tmax, transpose(m[1:S,1:tmax]))

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
fw = FoodWeb(nichemodel, 10, C = .33, check_cycle = true, check_disconnected = true)
max_trophic_level(fw.A)
A= fw.A

w = sim_int_mat(A;
            ρ = 0.0, alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            K = 10.0,
            σₑ = .6,
            Z = 10, h = 2.0, c = 0.0,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 5000, last = 500, dt = .1, return_sol = false, digits = 5);
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

