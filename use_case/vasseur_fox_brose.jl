println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using BEFWM2, LinearAlgebra, DifferentialEquations, DataFrames, CSV, Distributed, Distributions, ProgressMeter, BenchmarkTools

ncpu = length(Sys.cpu_info())
println("Using $(ncpu -2) cores")

#Flag enables all the workers to start on the project of the current dir
flag = "--project=~/xStressorsStabBEFW/"
println("Workers run with flag: $(flag)")
addprocs(ncpu - 2, exeflags=flag)

@everywhere import Pkg
@everywhere using DifferentialEquations, BEFWM2, Distributions, ProgressMeter
@everywhere include("../src/stochastic_mortality_model.jl")

import Random.seed!

seed!(22)

A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0]
foodweb = FoodWeb(A)

# Parameter product
#
# Along correlation between environmental stochasticity and strength of
# stochasticity
ρ = -1.0:0.1:1.0
σₑ = 0.0:.05:.6
rep = 1:5
Z_levels = 10 .^ [1:1:5;] # average predator-prey mass ratio
fr_types = [
    (h = 1.0, c = 0.0, type = "type II"),
    (h = 1.0, c = 1.0, type = "PI"),
    (h = 2.0, c = 0.0, type = "type III"),
]
names = (:ρ, :σₑ, :rep, :Z, :fr)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(ρ, σₑ, rep, Z_levels, fr_types))[:]


# Build the covariance matrix of the size of two times the number of species
S = BEFWM2.richness(foodweb)
vc = zeros(S * 2, S * 2)
# Only the consumer have a variance
# vc[diagind(vc)][S+2:S+3] .= σₑ
vc[diagind(vc)] .= 1.0
vc



@everywhere function mysim(A, n, Z, f, ρ, σₑ; max = 50000, last = 25000, dt = 0.1, corr_mat = vc, return_sol = false)

    fw = FoodWeb(A, Z = Z)
    S = BEFWM2.richness(fw)
    # change the parameters of the functional response
    # (we want the original functional response, as defined in Yodzis and Ines original paper)
    funcrep = BioenergeticResponse(fw; h = f.h, c = f.c)
    # generate the model parameters
    p = ModelParameters(fw; functional_response = funcrep, env_stoch = EnvStoch(σₑ))
    stoch_starting_val = [0; 0; 0; 0]
    u0 = [rand(size(A, 1)); stoch_starting_val]

    # Correlation between consumer
    corr_mat[S+2, S+3] = ρ
    corr_mat[S+3, S+2] = ρ

    # Generate the Wiener Process
    wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

    prob = SDEProblem(
                      mydBdt!,
                      stochastic_process,
                      u0,
                      [0, max],
                      p,
                      noise = wiener
                     )
    # Simulate
    timing = @elapsed m = try
        solve(prob;
              saveat = collect(0:1:max),
              dt = dt,
              adaptive = false
             )

    catch
        (t = 0, x = missing)
    end

    if return_sol
        return m
    end

    if length(m.t) == max + 1

        cv = foodweb_cv(m, last = last, idxs = [1, 2, 3, 4])
        sync_cons = foodweb_cv(m, last = last, idxs = [2, 3]).synchrony

        out = (
               rep = n,
               Z = Z,
               fr = f.type,
               ρ = ρ,
               σₑ = σₑ,
               sim_timing = timing,
               stab_com = 1 / cv.cv_com,
               avg_cv_sp = cv.avg_cv_sp,
               sync = cv.synchrony,
               stab_cons1 = 1 / cv.cv_sp[2],
               stab_cons2 = 1 / cv.cv_sp[3],
               stab_pred = 1 / cv.cv_sp[4],
               sync_cons = sync_cons
              )
    else
        out = (
               rep = n,
               Z = Z,
               fr = f.type,
               ρ = ρ,
               σₑ = σₑ,
               sim_timing = timing,
               stab_com = missing,
               avg_cv_sp = missing,
               sync = missing,
               stab_cons1 = missing,
               stab_cons2 = missing,
               stab_pred = missing,
               sync_cons = missing
              )
    end
    out
end

println("Test with batch size = 5")
@elapsed testres = @showprogress pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 1000, dt = .1, corr_mat = vc), param[10000:10010], batch_size = 5)
println("Test with batch size = 1")
@elapsed testres2 =  @showprogress pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 100, dt = .1, corr_mat = vc), param[10000:10010], batch_size = 1)

println("Real simulation")
@elapsed res = pmap(x -> mysim(A, x.rep, x.Z, x.fr, x.ρ, x.σₑ, max = 5000, last = 1000, dt = .1, corr_mat = vc), param, batch_size = 100)

df = DataFrame(res)

CSV.write("vasseur_fox__brose_res.csv", df)
