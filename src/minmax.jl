function local_extrema(sol; sp::Int = 1, tmin::Int = 0, tmax = nothing, nt::Int = 100000, verbose::Bool = true)
    if tmax == nothing
        tmax = length(sol.t)
    end

    # Interpolation of solution
    interp_sol = sol(range(tmin, tmax, length = nt), idxs = sp)
    bm = interp_sol.u

    l = []
    for i in 1:size(bm, 1)-2

        if i == 1
            d = bm[i+1] - bm[i]
            d1 = d 
        else
            d = bm[i+1] - bm[i]
            d1 = bm[i+2] - bm[i+1]
        end

        if sign(d1) != sign(d)
            t_ori = tmin + (1 + i) * (tmax - tmin) / nt
            push!(l, (t = interp_sol.t[i+1], extrema = bm[i+1]))
            if verbose == true
                println("t = $(interp_sol.t[i+1]), ", "extrema = $(bm[i + 1])")
            end
        end

    end
    out = (t = [l[i].t for i in 1:size(l, 1)],
           extrema = [l[i].extrema for i in 1:size(l, 1)]
          )
    return out 
end
