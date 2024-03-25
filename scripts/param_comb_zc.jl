using CSV, StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics, Random

sigma = .1:.1:0.6#[0.1, 0.3, 0.6]#.1:.2:0.6
Z = [1, 5, 10, 25, 50, 100]
ρ = 0:.25:1
c = [0.0, 1.0]

fw_d2 = DataFrame(Arrow.Table("scripts/fw_comb_ct_S_h_d2.arrow"))

rep = unique(fw_d2[:, :rep])
ct = unique(fw_d2[:, :C])

Random.seed!(1234)
rep_to_keep = sample(fw_d2[:, :rep], 15)
in_rep = in(rep_to_keep)
mask = in_rep.(fw_d2.rep)

df_fw = fw_d2[mask,:]


names = (:fw_id, :Z, :sigma, :rho, :c)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(df_fw[:, :fw_id], Z, sigma, ρ, c)
           )[:]

param_complete = innerjoin(DataFrame(param), df_fw, on = :fw_id)

# Cleaning
# SparseArrays to dense
param_complete[!, :A] = param_complete[:, :fw]
select!(param_complete, Not([:C, :fw]))
param_complete[!, :sim_id] = 1:nrow(param_complete)

nrow(param_complete) ÷ 2000

Arrow.write("scripts/param_comb_zc.arrow", param_complete)
ti = DataFrame(Arrow.Table("scripts/param_comb_zc.arrow"))

