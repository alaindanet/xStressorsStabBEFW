function sim_int_mat_check_disconnected(A;
        d = nothing,
        da = (ap = .4, ai = .4, ae = .4),
        ρ = 0.0, σₑ = .5,
        alpha_ij = 0.5, Z = 100,
        h = 2.0, c = 0.0, K = 10.0,
        r = 1.0,
        max = 5000, last = 1000,
        dt = 0.1,
        K_alpha_corrected = true,
        dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
        extinction_threshold = 1e-6,
        gc_thre = .02,
        return_sol = false,
        re_run = true,
        remove_disconnected = true,
        dt_rescue = 0.05,
        digits = nothing
    )

    ti = false
    i = 1
    starting_bm = nothing
    p = nothing
    println("Re-run until no more disconnected species \
            over the last $(last) timesteps.")
    while ti != true
         println("Iteration = $i.")
        if i == 1
            max = max
        else
            max = max
        end
        output = sim_int_mat(A;
                        B0 = starting_bm,
                        p  = nothing,
                        d = d,
                        da = da,
                        ρ = ρ, σₑ = σₑ,
                        alpha_ij = alpha_ij, Z = Z,
                        h = h, c = c, K = K,
                        r = r,
                        max = max, last = last,
                        dt = dt,
                        K_alpha_corrected = K_alpha_corrected,
                        dbdt = dbdt,
                        extinction_threshold = extinction_threshold,
                        gc_thre = gc_thre,
                        return_sol = false,
                        return_param = true,
                        re_run = re_run,
                        dt_rescue = dt_rescue,
                        digits = digits)
        if ismissing(output.omega) || ismissing(output.alive_species)
            println("No more species or problem in starting biomass.")
            return Base.structdiff(output, NamedTuple{(:param,)})
            break
        end
        # Retrieve network
        old_A = output.omega .> 0
        alive_species = output.alive_species
        # Check
        disconnected_species = check_disconnected_species(old_A, alive_species)
        ti = length(disconnected_species) == 0
        if ti == true
            println("No disconnected species, all fine.")
            # Remove model parameters from output
            return Base.structdiff(output, NamedTuple{(:param,)})
        end

        # Build the vector of biomasses
        biomass_vector = zeros(dim(old_A))
        biomass_vector[alive_species] = output.bm_sp
        killed_species = kill_disconnected_species(old_A;
                                                    alive_species = alive_species,
                                                    bm = biomass_vector)
        mask_sp_to_keep = killed_species .!= 0.0
        idxs = 1:dim(old_A)
        idxs_to_keep = idxs[mask_sp_to_keep]
        # Rebuilding network or set disconnected species to 0
        if remove_disconnected == true
            A = old_A[mask_sp_to_keep, mask_sp_to_keep]
            starting_bm = biomass_vector[mask_sp_to_keep]
            println("Rebuilding model without disconnected species.")
        else
            #starting_bm = biomass_vector[mask_sp_to_keep]
            #p = filter_model_parameters(output.param, idxs = idxs_to_keep)
            #A = p.network.A
            starting_bm = killed_species
            A = A
            println("Keep the same model without disconnected species.")
        end


        if length(starting_bm) != dim(A)
            println("Length of starting biomass differs from A: \
                    Starting biomass: $starting_bm, A: $(A), dim A: $(dim(A)).")
            println("alive species: $alive_species, killed species: $killed_species.")
        end
        i = i + 1
    end
    Base.structdiff(output, NamedTuple{(:param,)})
end

"""

# Examples/test

## Basic example 
ti = sim_int_mat([0 0; 0 0];
            ρ = 0,
            B0 = nothing,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = false,
            re_run = false,
            digits = 5
           )
# Wrong biomass length: 
ti = sim_int_mat([0 0; 0 0];
            ρ = 0,
            B0 = rand(1),
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = true,
            re_run = true,
            digits = 5
           )
# Empty biomass or A: 
empty_bm = rand(1)[[false]]
ti = sim_int_mat([0 0; 0 0];
            ρ = 0,
            B0 = empty_bm,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = true,
            re_run = true,
            digits = 5
           )
empty_A = [0 0; 0 0][[false, false], [false, false]]
ti = sim_int_mat(empty_A;
            ρ = 0,
            B0 = empty_bm,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = true,
            re_run = true,
            digits = 5
           )
ti = sim_int_mat(empty_A;
            ρ = 0,
            B0 = rand(2),
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = true,
            re_run = true,
            digits = 5
           )
# Negative bm
negative_bm = [0, -3]
ti = sim_int_mat([0 0; 0 0];
            ρ = 0,
            B0 = negative_bm,
            alpha_ij = 0.5,
            d = nothing,
            da = (ap = .4, ai = .4, ae = .4),
            σₑ = .1, Z = 10,
            h = 1.0, c = 0.0, K = 10,
            dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
            max = 1000, last = 100,
            K_alpha_corrected = true,
            dt = 0.1, gc_thre = .1,
            dt_rescue = .05,
            return_sol = false,
            re_run = true,
            digits = 5
           )
# NaN biomass
ti = sim_int_mat([0 0 0; 0 0 0; 1 1 0];
                 B0 = [NaN, NaN, NaN],
                 ρ = 1.0, alpha_ij = 0,
                 d = 0.0,
                 σₑ = .5, Z = 1000, h = 1.25, c = 0.0, K = 1.0,
                 dbdt = EcologicalNetworksDynamics.stoch_d_dBdt!,
                 max = 1000, last = 500, dt = 0.1, return_sol = true)
"""
function sim_int_mat(A;
        d = 0, ρ = 0.0, σₑ = .5,
        alpha_ij = 0, Z = 100,
        h = 2.0, c = 1.0, K = 1.0,
        r = 1.0,
        da = nothing,
        max = 5000, last = 1000,
        dt = 0.1,
        K_alpha_corrected = true,
        B0 = nothing,
        p = nothing,
        dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
        extinction_threshold = 1e-5,
        gc_thre = .02,
        return_sol = false,
        return_param = false,
        re_run = false,
        dt_rescue = 0.05,
        digits = nothing)

    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        GC.safepoint()
    end

    if !isnothing(B0)
        B0 = sanatize_biomass(B0)
    end

    if dim(A) != 0 && check_starting_bm(A, B0)
        if isnothing(p)

            fw = FoodWeb(A, Z = Z, quiet = true)
            S = richness(fw)
            C = connectance(fw)

            # Parameters of the functional response
            funcrep = BioenergeticResponse(fw; h = h, c = c)
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
        end

        ω = p.functional_response.ω
        if isnothing(B0)
            B0 = rand(S)
        end
        # Simulate
        timing = @elapsed m = try
            simulate(p, B0;
                     rho = ρ,
                     dt = dt,
                     tmax = max,
                     extinction_threshold = extinction_threshold,
                     diff_code_data = (dbdt, p),
                     verbose = false
                    );

        catch
            (t = 0, x = missing)
        end

        # If simulation did not fail
        if length(m.t) != 0
            # If unstable, it means that dt was too big
            if m.retcode == DifferentialEquations.ReturnCode.Unstable && !isnothing(dt_rescue)
                dt = dt_rescue
                println("rerun with dt = $dt_rescue.")
                # Simulate
                timing = @elapsed m = try
                    simulate(p, B0;
                             rho = ρ,
                             dt = dt,
                             tmax = max,
                             extinction_threshold = extinction_threshold,
                             diff_code_data = (dbdt, p),
                             verbose = false
                            );

                catch
                    (t = 0, x = missing)
                end
            end
        end


        # Re-run simulations until no more extinction
        if re_run
            if length(m.t) >= last
                ti = EcologicalNetworksDynamics.check_last_extinction(m;
                                                                      idxs = 1:S,
                                                                      last = last)
                i = 1
                while ti != true
                    B0 = m[1:S,end]
                    println("Re-run until no more extinctions \
                            over the last $(last) timesteps. Iteration = $i.")

                    #u0 = m[S+1:S*2,end]
                    p = get_parameters(m)
                    m = simulate(p, B0;
                                 rho = ρ,
                                 dt = dt,
                                 tmax = round(Int, last + 500),
                                 extinction_threshold = extinction_threshold,
                                 diff_code_data = (dbdt, p),
                                 verbose = false
                                );
                    ti = EcologicalNetworksDynamics.check_last_extinction(m;
                                                                          idxs = 1:S,
                                                                          last = last)
                    if ti == true
                        break
                    end
                    i = i + 1
                end
            end
        end

        if return_sol
            return m
        end
    else
        println("No more species, returns missing.")
        m = (t = 0, x = missing)
        ω = missing
        d = missing
        S = 0
        C = 0
        timing = missing
    end

    out = get_stab_fw(m; last = last, digits = digits)
    # Collect timeseries, only alive species
    time_series = get_time_series(m; last = last,
                                  idxs = out.alive_species,
                                  digits = digits)
    if length(m.t) >= last
        ω = get_parameters(m).functional_response.ω
        #ω = ω[out.alive_species, out.alive_species]
        d = get_parameters(m).biorates.d
        d = d[out.alive_species]
        if !isnothing(digits)
            ω = round.(ω, digits = digits)
            d = round.(d, digits = digits)
        end
    end

    out = merge(
                (S = S, ct = C, rho = ρ, env_stoch = σₑ,
                 Z = Z, omega = ω, d = d, timing = timing),
                out,
                time_series
               )
    if return_param
        param = try
            get_parameters(m)
        catch
            missing
        end
        out = merge(out, (param = param,))
    end

    out
end

function get_time_series(m; last = 10, idxs = nothing, digits = 5)
    names_output = (
                    :species,
                    :stoch
                   )
    if length(m.t) >= last
        S = richness(get_parameters(m).network)
        if isnothing(idxs)
            idxs = 1:S
        end
        values = (
                  transpose(m[1:S, end-(last-1):end]),
                  transpose(m[S+1:2*S, end-(last-1):end])
                 )
        values = map(x -> x[:, idxs], values)
        if !isnothing(digits)
            values = map(x -> round.(x, digits = digits), values)
        end
    else
        values = repeat([missing], length(names_output))
    end
    (; zip(names_output, values)...)
end

function get_stab_fw(m; last = 10, digits = nothing, kwargs...)
    names_output = (
                    :richness          ,
                    :stab_com          ,
                    :avg_cv_sp         ,
                    :sync              ,
                    :total_biomass     ,
                    :bm_sp             ,
                    :cv_sp             ,
                    :alive_species     ,
                    :tlvl              ,
                    :w_avg_tlvl        ,
                    :max_tlvl          ,
                    :int_strength      ,
                    :avg_int_strength  ,
                    :max_int_strength  ,
                    :min_int_strength  ,
                    :gini_int_strength ,
                    :max_int           ,
                    :avg_max_int       ,
                    :max_max_int       ,
                    :min_max_int       ,
                    :gini_max_int      ,
                    :omnivory          ,
                    :avg_omnivory
                   )

    if length(m.t) >= last && maximum(m.t) != 0
        p = get_parameters(m)
        fw = p.network
        bm_cv = coefficient_of_variation(m; last = last, kwargs...)
        bm = biomass(m, last = last)
        emp_int_strength = empirical_interaction_strength(m, p, last = last)
        non_zero_int = [i for i in emp_int_strength.mean if i != 0]
        max_int = max_interaction_strength(p)
        non_zero_max_int = [i for i in max_int if i != 0]
        troph_struc = trophic_structure(m, last = last)
        alive_species = troph_struc.alive_species
        omnivory_alive = omnivory(emp_int_strength.mean[troph_struc.alive_species, troph_struc.alive_species])

        values = [
                  richness(m, last = last),
                  1 / bm_cv.community,
                  bm_cv.average_species,
                  bm_cv.synchrony,
                  bm.total,
                  bm.species[alive_species],
                  bm_cv.species,
                  alive_species,
                  troph_struc.alive_trophic_level,
                  troph_struc.weighted_average,
                  troph_struc.max,
                  emp_int_strength.mean,
                  mean(non_zero_int),
                  maximum(non_zero_int, init = 0),
                  minimum(non_zero_int, init = 0),
                  gini(non_zero_int),
                  max_int,
                  mean(non_zero_max_int),
                  maximum(non_zero_max_int, init = 0),
                  minimum(non_zero_max_int, init = 0),
                  gini(non_zero_max_int),
                  omnivory_alive,
                  mean(omnivory_alive)
                 ]
        if !isnothing(digits)
            alive_species_mask = [names_output[i] == :alive_species for i in 1:length(names_output)]
            values = map(x -> round.(x, digits = digits), values)
            values[findall(alive_species_mask)[1]] = round.(Int, values[findall(alive_species_mask)[1]])
        end
    else
        values = repeat([missing], length(names_output))
    end
    (; zip(names_output, values)...)
end

function mysim(A, n, Z, f, ρ, σₑ; max = 50000, last = 25000, dt = 0.1, corr_mat = vc, return_sol = false)

    fw = FoodWeb(A, Z = Z)
    S = richness(fw)
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

        cv = cv(m, last = last, idxs = [1, 2, 3, 4])
        sync_cons = cv(m, last = last, idxs = [2, 3]).synchrony

        out = (
               rep        = n,
               Z          = Z,
               fr         = f.type,
               ρ          = ρ,
               σₑ         = σₑ,
               sim_timing = timing,
               stab_com   = 1 / cv.cv_com,
               avg_cv_sp  = cv.avg_cv_sp,
               sync       = cv.synchrony,
               stab_res   = 1 / cv.cv_sp[1],
               stab_cons1 = 1 / cv.cv_sp[2],
               stab_cons2 = 1 / cv.cv_sp[3],
               stab_pred  = 1 / cv.cv_sp[4],
               sync_cons  = sync_cons,
               bm_res     = cv.bm_sp[1],
               bm_cons1   = cv.bm_sp[2],
               bm_cons2   = cv.bm_sp[3],
               bm_pred    = cv.bm_sp[4]
              )
    else
        out = (
               rep        = n,
               Z          = Z,
               fr         = f.type,
               ρ          = ρ,
               σₑ         = σₑ,
               sim_timing = timing,
               stab_com   = missing,
               avg_cv_sp  = missing,
               sync       = missing,
               stab_res   = missing,
               stab_cons1 = missing,
               stab_cons2 = missing,
               stab_pred  = missing,
               sync_cons  = missing,
               bm_res     = missing,
               bm_cons1   = missing,
               bm_cons2   = missing,
               bm_pred    = missing
              )
    end
    out
end

"""
kill_disconnected_species(A; alive_species = nothing, bm = nothing)

Returns a vector of species biomass, all the disconnected and extinct species being set to a
biomass of 0.

# Arguments:

  - `A`: Adjacency matrix
  - `alive_species`: vector of alive species indices. See [`living_species`](@ref)
  - `bm`: vector of species biomass

# Examples

ti = [
 0 0 0 0;
 0 0 0 0;
 1 0 0 0;
 0 1 0 0
]

bm = [1, 1, 3, 4]
alive = [1, 3, 4]
kill_disconnected_species(ti, alive_species = [1, 2, 3, 4], bm = bm) == [1, 1, 3, 4]
kill_disconnected_species(ti, alive_species = [], bm = bm) == [0, 0, 0, 0]
kill_disconnected_species(ti, alive_species = [1, 3, 4], bm = bm) == [1, 0, 3, 0]
kill_disconnected_species(ti, alive_species = [1, 2, 4], bm = bm) == [0, 1, 0, 4]
kill_disconnected_species(ti, alive_species = [1, 2, 4], bm = bm) == [0, 1, 0, 4]

# Filter species biomass and then adjacency matrix

## Case 1
to = kill_disconnected_species(ti, alive_species = [1, 2, 3, 4], bm = bm)
mask = to .== 0
A = ti[mask, mask]
new_bm = bm[mask]
dim(A) == length(new_bm)

## Case 2
to = kill_disconnected_species(ti, alive_species = [1, 2, 4], bm = bm)
mask = to .== 0
A = ti[mask, mask]
new_bm = bm[mask]
dim(A) == length(new_bm)


"""
function kill_disconnected_species(A; alive_species = nothing, bm = nothing)
    bm = deepcopy(bm)
    disconnected_alive_species = check_disconnected_species(A, alive_species)
    in_disconnected = in(disconnected_alive_species)
    mask_alive_disconnected = in_disconnected.(alive_species)
    alive_to_keep = alive_species[.!mask_alive_disconnected]
    bm_to_set_to_zero = 1:length(bm) .∉ [alive_to_keep]
    bm[bm_to_set_to_zero] .= 0
    bm
end

"""
    remove_disconnected_species(A, alive_species)

# Examples

ti = [
 0 0 0 0;
 0 0 0 0;
 1 0 0 0;
 0 1 0 0
]
remove_disconnected_species(ti, [1, 2, 3])
"""
function remove_disconnected_species(A, alive_species)

    disconnected_alive_species = check_disconnected_species(A, alive_species)
    in_disconnected = in(disconnected_alive_species)
    mask_alive_disconnected = in_disconnected.(alive_species)
    to_keep = .!(mask_alive_disconnected)

    living_A = A[alive_species, alive_species]
    living_A[to_keep, to_keep]

end

"""
    check_disconnected_species(A, alive_species)

# Examples

ti = [
 0 0 0 0;
 0 0 0 0;
 1 0 0 0;
 0 1 0 0
]

check_disconnected_species(ti, [1, 2, 3]) == [2]
check_disconnected_species(ti, [1, 2, 4]) == [1]
check_disconnected_species(ti, [2, 3, 4]) == [3]
check_disconnected_species(ti, [1, 3, 4]) == [4]
check_disconnected_species(ti, [1, 2, 3, 4]) == []

"""
function check_disconnected_species(A, alive_species; verbose = false)
    # Get binary matrix
    A = A .> 0
    living_A = A[alive_species, alive_species]
    cons = sum.(eachrow(A)) .!= 0
    prod = sum.(eachrow(A)) .== 0


    alive_cons = cons[alive_species]
    alive_prod = prod[alive_species]

    species_with_no_pred = sum.(eachcol(living_A)) .== 0
    species_with_no_prey = sum.(eachrow(living_A)) .== 0

    disconnected_prod = alive_prod .&& species_with_no_pred
    disconnected_cons = alive_cons .&& species_with_no_prey

    if sum([disconnected_prod; disconnected_cons]) > 0 & verbose
        println("There are $(sum(disconnected_prod)) disconnected producers
            and $(sum(disconnected_cons)) consumers.")
    end

    alive_species[ disconnected_prod .|| disconnected_cons ]
end


function check_starting_bm(A, x)
    out = isnothing(x) || (!isnothing(x) && length(x) != 0 && length(x) == dim(A))
    if !isnothing(x) && length(x) != 0 && length(x) != dim(A)
        println("Starting biomass vector length do not fit species number in the network.\
                dim(A): $(dim(A)), length(biomass) = $(length(x)).
                ")
    end
    out
end

function sanatize_biomass(bm)
    # For negative biomass
    mask = bm .< 0
    if any(mask)
        println("Some starting biomass are negative: $(bm[mask]), put them to 0.")
    end
    bm[mask] .= 0
    bm
end

function filter_model_parameters(p; idxs = nothing)
    # Extract body mass, network, metabolic_class
    fw_new = FoodWeb(p.network.A[idxs, idxs],
                     M = p.network.M[idxs],
                     metabolic_class = p.network.metabolic_class[idxs]
                    )

    # BioEnergeticResponse
    # Extract omega and K, keep same alpha
    new_bioner = BioenergeticResponse(fw_new,
                                      h = p.functional_response.h,
                                      c = p.functional_response.c[idxs],
                                      B0 = p.functional_response.B0[idxs],
                                      ω = p.functional_response.ω[idxs, idxs]
                                     )
    # Environment
    new_env = Environment(fw_new, K = p.environment.K[idxs])

    # BioRates 
    new_biorates = BioRates(fw_new;
                            r = p.biorates.r[idxs],
                            d = p.biorates.d[idxs]
            )
    # ProducerCompetition
    new_prod = ProducerCompetition(fw_new, α = p.producer_competition.α[idxs, idxs])

    # EnvStoch
    env_stoch = p.env_stoch

    ModelParameters(
                    fw_new,
                    biorates = new_biorates,
                    functional_response = new_bioner,
                    environment = new_env,
                    producer_competition = new_prod,
                    env_stoch = env_stoch
                   )
end

function simCS(C, S;
        Z = 100,
        d = 0, σₑ = .5, ρ = 0.0,
        h = 2.0, c = 0.0,
        r = 1.0, K = 5.0, alpha_ij = 0.5,
        max = 5000, last = 1000, dt = 0.1,
        return_sol = false,
        dbdt = EcologicalNetworksDynamics.stoch_m_dBdt!,
        K_alpha_corrected = true,
        Ctol = .02,
        gc_thre = .02
    )

    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        GC.safepoint()
    end

    # Generate a food-web with a given connectance
    fw = try
        FoodWeb(nichemodel, S, C = C, Z = Z, tol_C = Ctol,
                check_cycle = true,
                check_disconnected = true)
    catch
        missing
    end

    if ismissing(fw)
        println("FoodWeb generation has failed with C = $(C), S = $(S)")
        A = missing
    else
        A = fw.A
    end

    # If the generation of the food-web worked
    if !ismissing(fw)
        # Parameters of the functional response
        funcrep = BioenergeticResponse(fw; h = h, c = c)
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
        # generate the model parameters
        p = ModelParameters(fw;
                            biorates = BioRates(fw; r = r, d = d),
                            functional_response = funcrep,
                            environment = env,
                            producer_competition = ProducerCompetition(fw, αii = 1.0, αij = alpha_ij),
                            env_stoch = EnvStoch(σₑ))
        ω = p.functional_response.ω
    B0 = rand(S)
    # Simulate
    timing = @elapsed m = try
        simulate(p, B0;
                 rho = ρ,
                 dt = dt,
                 tmax = max,
                 extinction_threshold = 1e-5,
                 diff_code_data = (dbdt, p),
                 verbose = false
                );
        catch
            (t = 0, x = missing)
        end

    else
        m = (t = 0, x = missing)
        timing = missing
        ω = missing
    end

    if return_sol
        return m
    end

    out = get_stab_fw(m; last = last)
    # Collect timeseries
    time_series = get_time_series(m; last = last)

    out = merge(
                (S = S, ct = C, omega = ω, rho = ρ, env_stoch = σₑ, Z = Z, timing = timing),
                out,
                time_series
               )
    out
end

