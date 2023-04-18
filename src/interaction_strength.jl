"""
    max_interaction_strength(p::ModelParameters;)

Returns the theoretical maximum interaction strengths from consumers to resources, realized
when the density of a given resource approaches 0.

# References

McCann, K., Hastings, A., & Huxel, G. R. (1998). Weak trophic interactions and the balance
of nature. Nature, 395(6704), Art. 6704. https://doi.org/10.1038/27427
"""
function max_interaction_strength(ω, x, y, B0, h;)
    S = size(ω, 2)
    int = zeros(S, S)

    for i in 1:S
        int[i, :] = [(ω[i, j] * x[i] * y[i]) / ((B0[i])^h) for j in 1:S]
    end

    int
end

function max_interaction_strength(p::ModelParameters;)
    max_interaction_strength(p.functional_response.ω,
                             p.biorates.x,
                             p.biorates.y,
                             p.functional_response.B0,
                             p.functional_response.h
                            )
end

"""
    gini(x)

Computes the Gini coefficient, a index of inequality. A Gini coefficient of 1 indicates
complete inequality while a value of 0 indicates complete equality.

# References

https://gist.github.com/ZacLN/fd34052b6a17ea9274f1e92cf82d6472
https://en.wikipedia.org/wiki/Gini_coefficient

"""
function gini(x)
    n = length(x)
    xx = sort(x)
    2*(sum(collect(1:n).*xx))/(n*sum(xx))-1
end

"""
    empirical_interaction_strength(solution, params; kwargs...)

Computes the empirical trophic interaction strength over the `last` timesteps.
It returns the mean, max, min, standard deviation of each trophic interaction in the
network, and return the timeseries as well.

"""
function empirical_interaction_strength(solution, params::ModelParameters; kwargs...)

    measure_on = extract_last_timesteps(solution; kwargs...)

    S = richness(params.network)
    ntimestep = size(measure_on, 2)
    out = zeros(S, S, ntimestep)
    for i in 1:ntimestep
        out[:, :, i] = empirical_interaction_strength(measure_on[:, i], params)
    end

    (
     mean = mean(out, dims = 3)[:, :, 1],
     max = maximum(out, dims = 3)[:, :, 1],
     min = minimum(out, dims = 3)[:, :, 1],
     std = std(out, dims = 3)[:,:, 1],
     all = out
   )

end

function empirical_interaction_strength(B::Vector{Float64}, params::ModelParameters)

    S = size(params.network.species, 1)
    int = zeros(S, S)
    h = params.functional_response.h
    B0 = params.functional_response.B0
    ω = params.functional_response.ω
    x = params.biorates.x
    y = params.biorates.y
    c = params.functional_response.c

    for i in 1:S
        int[i, :] = [
                     ( x[i] * y[i] * B[i] * ω[i, j] * (B[j])^h ) /
                     ( (B0[i])^h + c[i] * B[i] * (B0[i])^h  + sum(ω[i,:] .* (B .^h)))
                     for j in 1:S
                    ]
    end
    sparse(int)
end
