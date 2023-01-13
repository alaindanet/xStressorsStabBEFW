function simCS(C, S, Z, h, c, σₑ, K; max = 50000, last = 25000, dt = 0.1, return_sol = false)

    # Generate a food-web with a given connectance
    fw = try
        global tmp_fw
        ct = 0
        while ct != C
            tmp_fw = FoodWeb(nichemodel, S, C = C, Z = Z)
            ct = round(connectance(tmp_fw), digits = 2)
        end
        tmp_fw
    catch
        missing
    end

    if ismissing(fw)
        println("FoodWeb generation has failed with C = $(C), S = $(S)")
    end

    # If the generation of the food-web worked
    if !ismissing(fw)
        # Parameters of the functional response
        funcrep = BioenergeticResponse(fw; h = h, c = c)
        # generate the model parameters
        p = ModelParameters(fw;
                            functional_response = funcrep,
                            environment = Environment(fw, K = K/length(BEFWM2.producers(fw))),
                            env_stoch = EnvStoch(σₑ))
        stoch_starting_val = repeat([0], S)
        u0 = [rand(S); stoch_starting_val]

        # Make the stochastic matrix
        corr_mat = Diagonal(repeat([1.0], S * 2))
        # Generate the Wiener Process
        wiener = CorrelatedWienerProcess(corr_mat, 0.0, zeros(size(corr_mat, 1)))

        prob = SDEProblem(
                          mydBdt!,
                          gen_stochastic_process,
                          u0,
                          [0, max],
                          p,
                          noise = wiener
                         )
        # Try to simulate
        timing = @elapsed m = try
            solve(prob;
                  saveat = collect(0:1:max),
                  dt = dt,
                  adaptive = false
                 )

        catch
            (t = 0, x = missing)
        end

    else
        m = (t = 0, x = missing)
        timing = missing
    end

    if return_sol
        return m
    end

    if length(m.t) == max + 1
        bm_cv = cv(m, last = last, idxs = collect(1:1:S))
        bm = biomass(m, last = last, idxs = collect(1:1:S))
        int_strength = empirical_interaction_strength(m, p, last = last)
        non_zero_int = [i for i in int_strength if i != 0]
        max_int = max_interaction_strength(p)
        non_zero_max_int = [i for i in max_int if i != 0]
        tlvl = trophic_levels(fw)

        out = (
               richness          = S,
               ct                = C,
               Z                 = Z,
               env_stoch         = σₑ,
               sim_timing        = timing,
               stab_com          = 1 / bm_cv.cv_com,
               avg_cv_sp         = bm_cv.avg_cv_sp,
               sync              = bm_cv.synchrony,
               total_biomass     = bm.total,
               bm_sp             = bm.sp,
               cv_sp             = bm_cv.cv_sp,
               tlvl              = tlvl,
               w_avg_tlvl        = sum(tlvl .* (bm.sp ./ sum(bm.sp))),
               max_tlvl          = maximum(tlvl),
               int_strength      = int_strength,
               avg_int_strength  = mean(non_zero_int),
               max_int_strength  = maximum(non_zero_int),
               min_int_strength  = minimum(non_zero_int),
               gini_int_strength = gini(non_zero_int),
               max_int           = max_int,
               avg_max_int       = mean(non_zero_max_int),
               max_max_int       = maximum(non_zero_max_int),
               min_max_int       = minimum(non_zero_max_int),
               gini_max_int      = gini(non_zero_max_int)
              )
    else
        out = (
               richness          = S,
               ct                = C,
               Z                 = Z,
               env_stoch         = σₑ,
               sim_timing        = timing,
               stab_com          = missing,
               avg_cv_sp         = missing,
               sync              = missing,
               total_biomass     = missing,
               bm_sp             = missing,
               cv_sp             = missing,
               tlvl              = missing,
               w_avg_tlvl        = missing,
               max_tlvl          = missing,
               int_strength      = missing,
               avg_int_strength  = missing,
               max_int_strength  = missing,
               min_int_strength  = missing,
               gini_int_strength = missing,
               max_int           = missing,
               avg_max_int       = missing,
               max_max_int       = missing,
               min_max_int       = missing,
               gini_max_int      = missing
              )
    end
    out
end

function mysim(A, n, Z, f, ρ, σₑ; max = 50000, last = 25000, dt = 0.1, corr_mat = vc, return_sol = false)

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

