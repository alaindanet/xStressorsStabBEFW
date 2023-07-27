using CSV, StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics


param = DataFrame(Arrow.Table("param_comb_ct_S_h.arrow"))

param.sigma = replace(param.sigma, 0.1 => 0.5)
param.sigma = replace(param.sigma, 0.3 => 1)
param.sigma = replace(param.sigma, 0.6 => 1.5)

unique(param.sigma)

param.sim_id = 1:nrow(param)
# sim_id as 1st column:
select!(param, append!(["sim_id"], names(param)[names(param) .!= "sim_id"]))

# Write data.frame
Arrow.write("param_comb_ct_S_h_d_allometric.arrow", param)
