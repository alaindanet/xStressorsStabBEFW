println("Rundir is $(pwd())")

import Pkg
using Distributed, Serialization

first_sim = parse(Int, ARGS[1])
last_sim = parse(Int, ARGS[2])

println("Running parameters combination from $(ARGS[1]) to $(ARGS[2])")

#ncpu = maximum([length(Sys.cpu_info()), 15])
ncpu = 15

#Dir
proj_dir = expanduser("~/xStressorsStabBEFW")

#Flag enables all the workers to start on the project of the current dir
flag = "--project=" * proj_dir
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -1) cores")


@everywhere import Pkg, Random.seed!
@everywhere using DifferentialEquations, EcologicalNetworksDynamics, SparseArrays
@everywhere using LinearAlgebra, DataFrames
@everywhere using Distributions, ProgressMeter
@everywhere using StatsBase, CSV, Arrow
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../src/interaction_strength.jl")


param = DataFrame(Arrow.Table(joinpath(proj_dir, "scripts/param_comb_ct_S_h.arrow")))

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
println("Running warmup:K = $(pm.K), σₑ = $(pm.sigma), ρ = $(pm.rho)")

warmup = sim_int_mat([0 0; 0 0];
            ρ = pm.rho, alpha_ij = 0,
            d = 0.0,
            σₑ = pm.sigma, Z = pm.Z, h = pm.h, c = 0.0, K = pm.K,
            dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
            max = 50, last = 10, dt = 0.1, return_sol = false)
println("$(warmup)")


if last_sim > size(param, 1)
    last_sim = size(param, 1)
end
println("Running param sim from lines $first_sim to $last_sim")

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (fw_id = p.fw_id, productivity = p.K, h = p.h),
                               sim_int_mat(p.A;
                                           ρ = p.rho,
                                           alpha_ij = 0.5,
                                           d = 0.0,
                                           σₑ = p.sigma, Z = p.Z,
                                           h = p.h, c = 0.0, K = p.K,
                                           dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
                                           max = 5000, last = 100,
                                           K_alpha_corrected = true,
                                           dt = 0.1, gc_thre = .1,
                                           return_sol = false
                                          )
                              ),
                         param[first_sim:last_sim],
                         batch_size = 100
                        )

df = DataFrame(sim)
println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

file_ts = string("simCSh_", first_sim, "_", last_sim, "_ts.arrow")
Arrow.write(file_ts, select(df, [:species, :stoch]))
file = string("simCSh_", first_sim, "_", last_sim, ".arrow")
Arrow.write(file, select(df, Not([:species, :stoch, :cv_sp])))
