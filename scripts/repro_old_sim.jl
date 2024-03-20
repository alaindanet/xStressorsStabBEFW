include("../scripts/sim_setup.jl")

# Prepare saving
dest_dir = "/mnt/parscratch/users/bi1ahd/sim/simCSh_allo_repro_old/"

if !isdir(dest_dir)
    mkdir(dest_dir)
end

param = DataFrame(Arrow.Table(joinpath(proj_dir, "scripts/param_comb_ct_S_h_d2.arrow")))

# Reshape interaction matrix
reshape_array(vec) = reshape(vec, (
                                   round(Int, sqrt(length(vec))),
                                   round(Int, sqrt(length(vec)))
                                  )
                            )
param[!, :A] = map(x -> reshape_array(x), param[:, :A])

# Make a tuple vector
param = NamedTuple.(eachrow(param))

# Warm-up
pm = sample(param)
println("Running warmup:σₑ = $(pm.sigma), ρ = $(pm.rho)")

warmup = sim_int_mat([0 0; 0 0];
            ρ = pm.rho, alpha_ij = 0,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = pm.sigma, Z = pm.Z, h = pm.h, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 50, last = 10, dt = 0.1, return_sol = false)
println("$(warmup)")


if last_sim > size(param, 1)
    last_sim = size(param, 1)
end
println("Running param sim from lines $first_sim to $last_sim")

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (sim_id = p.sim_id, fw_id = p.fw_id, h = p.h),
                               sim_int_mat_check_disconnected(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = nothing,
                                           da = (ap = .4, ai = .4, ae = .4),
                                           σₑ = p.sigma, Z = p.Z,
                                           h = 2.0, c = 0.0, K = 10,
                                           dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                                           max = 2000, last = 500,
                                           K_alpha_corrected = true,
                                           dt = 0.1, gc_thre = .1,
                                           dt_rescue = 0.05,
                                           extinction_threshold = 1e-6,
                                           return_sol = false,
                                           remove_disconnected = true,
                                           re_run = true, # works only if you get rid of disconnected species
                                           digits = 5
                                          )
                              ),
                         param[first_sim:last_sim],
                         batch_size = 100
                        )

df = DataFrame(sim)
println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

# Saving
file_ts = string("simCSh_allo_d", first_sim, "_", last_sim, "_ts.arrow")
Arrow.write(joinpath(dest_dir, file_ts), select(df, [:sim_id, :species, :stoch]))

file = string("simCSh_allo_d", first_sim, "_", last_sim, ".arrow")
Arrow.write(joinpath(dest_dir, file), select(df, Not([:species, :stoch, :cv_sp])))
