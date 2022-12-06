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

        cv = foodweb_cv(m, last = last, idxs = [1, 2, 3, 4])
        sync_cons = foodweb_cv(m, last = last, idxs = [2, 3]).synchrony

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

