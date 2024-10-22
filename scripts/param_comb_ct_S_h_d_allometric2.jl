using CSV, StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics

nrep = 5
sigma = .1:.1:0.6#[0.1, 0.3, 0.6]#.1:.2:0.6
Z = [1, 5, 10, 25, 50, 100]
ρ = 0:.25:1
h = 1.0:.25:2.0#[2.0, 3.0]
S = [10, 20, 30, 40, 50]
C = 0.02:.02:.32

# 
names = (:rep, :S, :C, :Z, :sigma, :rho, :h)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(1:nrep, S, C, Z, sigma, ρ, h)
           )[:]

########################
#  Generate Food-webs  #
########################

rep = 1:1
names = (:rep, :richness, :connectance)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(rep, S, C)
           )[:]

fw_test = map(p -> try
                  FoodWeb(nichemodel, p.richness, C = p.connectance,
                          Z = 100, tol_C = .05,
                          check_cycle = true,
                          check_disconnected = true)
              catch
                  missing
              end,
              param
             )


param[.! ismissing.(fw_test)]
fw_test[.! ismissing.(fw_test)]


fw_test2 = map(p -> try
                  FoodWeb(nichemodel, p.richness, C = p.connectance,
                          Z = 100, tol_C = .02,
                          check_cycle = true,
                          check_disconnected = true)
              catch
                  missing
              end,
              param[.! ismissing.(fw_test)]
             )

param2 = param[.! ismissing.(fw_test)][.! ismissing.(fw_test2)]

# Try to get same number of sampling per richness
df = groupby(DataFrame(param2), :richness)
rounded_mean(data_col) = round(mean(data_col), digits = 2)
df2 = combine(df, :connectance => minimum, :connectance => rounded_mean, :connectance => maximum)
df3 = stack(df2, 2:4, value_name= :connectance)
select!(df3, [:richness, :connectance])

function try_foodweb(S; C = .1, tol_C = .05, n = 5, kwargs...)
    fw = (A = missing,)
    local i = 1

    while all([ismissing(fw.A), i <= n])
        println("i = $i")
        fw = try FoodWeb(nichemodel,
                         S;
                         C = C,
                         tol_C = tol_C,
                         check_cycle = true,
                         check_disconnected = true,
                         kwargs...
                        )
        catch
            (A = missing,)
        end
        i = i + 1
    end
    fw
end

fw = try_foodweb(10; C = 0.5,
            tol_C = .05,
            check_cycle = true,
            check_disconnected = true).A


df4 = repeat(df3, nrep)
df4[!, :rep] = reduce(vcat, [repeat([i], nrow(df3)) for i in 1:nrep])

foodweb = map(p -> (rep = p.rep, S = p.richness, C = p.connectance,
                    fw = try_foodweb(p.richness; C = p.connectance,
                                tol_C = .05,
                                check_cycle = true,
                                check_disconnected = true).A
                  ),
             NamedTuple.(eachrow(df4))
            )
df_fw = DataFrame(foodweb)
df_fw = df_fw[.! ismissing.(df_fw.fw),:]
# Create a foodweb_id
df_fw[!, :fw_id] = 1:nrow(df_fw)
Arrow.write("../fw_comb_ct_S_h_d3.arrow", df_fw)

################################
#  Other parameter combination #
################################

names = (:fw_id, :Z, :sigma, :rho, :h)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(df_fw[:, :fw_id], Z, sigma, ρ, h)
           )[:]

param_complete = innerjoin(DataFrame(param), df_fw, on = :fw_id)

# Cleaning
# SparseArrays to dense
param_complete[!, :A] = map(x -> Matrix(x), param_complete[:, :fw])
select!(param_complete, Not([:C,:rep, :fw]))
param_complete[!, :sim_id] = 1:nrow(param_complete)


Arrow.write("../param_comb_ct_S_h_d3.arrow", param_complete)


# Change Z gradient
ti = DataFrame(Arrow.Table("scripts/param_comb_ct_S_h_d3.arrow"))

#Z = [1, 5, 10, 25, 50, 100]
ti.Z = replace(ti.Z, 5 => 1000, 25 => 10000, 50 => 500)
unique(ti.Z)
Arrow.write("scripts/param_comb_ct_S_h_d4.arrow", ti)
ti = DataFrame(Arrow.Table("scripts/param_comb_ct_S_h_d4.arrow"))

