function simCS(C, S;
        Z = 100,
        d = 0, σₑ = .5, ρ = 0.0,
        h = 2.0, c = 1.0,
        r = 1.0, K = 1.0, alpha_ij = 0.5,
        max = 5000, last = 1000, dt = 0.1,
        return_sol = false,
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
    B0 = rand(S)
    # Simulate
    timing = @elapsed m = try
        simulate(p, B0;
                 rho = ρ,
                 dt = dt,
                 tmax = max,
                 extinction_threshold = 1e-5,
                 verbose = false
                );
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

    out = get_stab_fw(m; last = last)
    # Collect timeseries
    time_series = get_time_series(m; last = last)

    out = merge(
                (S = S, ct = C, A = A, rho = ρ, env_stoch = σₑ, Z = Z, timing = timing),
                out,
                time_series
               )
    out
end

function sim_int_mat(A;
        d = 0, ρ = 0.0, σₑ = .5,
        alpha_ij = 0, Z = 100,
        h = 2.0, c = 1.0, K = 1.0,
        r = 1.0,
        max = 5000, last = 1000,dt = 0.1,
        K_alpha_corrected = true,
        return_sol = false)

    fw = FoodWeb(A, Z = Z, quiet = true)
    S = richness(fw)

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
                        env_stoch = EnvStoch(σₑ)
                       )
    B0 = rand(S)
    # Simulate
    timing = @elapsed m = try
        simulate(p, B0;
                 rho = ρ,
                 dt = dt,
                 tmax = max,
                 extinction_threshold = 1e-5,
                 verbose = false
                );

    catch
        (t = 0, x = missing)
    end

    if return_sol
        return m
    end

    out = get_stab_fw(m; last = last)
    # Collect timeseries
    time_series = get_time_series(m; last = last)

    out = merge(
                (S = S, rho = ρ, env_stoch = σₑ, Z = Z, timing = timing),
                out,
                time_series
               )
    out
end

function get_time_series(m; last = 10)
    names_output = (
                    :species,
                    :stoch
                   )
    if length(m.t) >= last
        S = richness(get_parameters(m).network)
        values = (
                  transpose(round.(m[1:S, end-(last-1):end], digits = 5)),
                  transpose(round.(m[S+1:2*S, end-(last-1):end], digits = 5))
                 )
    else
        values = repeat([missing], length(names_output))
    end
    (; zip(names_output, values)...)
end

function get_stab_fw(m; last = 10)
    names_output = (
                    :richness          ,
                    :stab_com          ,
                    :avg_cv_sp         ,
                    :sync              ,
                    :total_biomass     ,
                    :bm_sp             ,
                    :cv_sp             ,
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
                    :gini_max_int,
                    :omnivory
                   )

    if length(m.t) >= last
        p = get_parameters(m)
        fw = p.network
        bm_cv = coefficient_of_variation(m, last = last)
        bm = biomass(m, last = last)
        emp_int_strength = empirical_interaction_strength(m, p, last = last)
        non_zero_int = [i for i in emp_int_strength.mean if i != 0]
        max_int = max_interaction_strength(p)
        non_zero_max_int = [i for i in max_int if i != 0]
        troph_struc = trophic_structure(m, last = last)

        values = (
                  richness(m, last = last),
                  1 / bm_cv.community,
                  bm_cv.average_species,
                  bm_cv.synchrony,
                  bm.total,
                  bm.species,
                  bm_cv.species,
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
                  omnivory(p; weighted = true)
                 )
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

