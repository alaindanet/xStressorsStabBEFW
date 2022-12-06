function max_interaction_strength(ω, x, y, B0, h;)
    S = size(ω, 2)
    int = zeros(S, S)

    for i in 1:S
        int[i, :] = [(ω[i, j] * x[i] * y[i]) / ((B0[i])^h) for j in 1:S]
    end

    int
end

function max_interaction_strength(p::ModelParameters;)
    max_interaction_strength(p.functional_response.ω, p.biorates.x, p.biorates.y, p.functional_response.B0, p.functional_response.h)
end

# https://gist.github.com/ZacLN/fd34052b6a17ea9274f1e92cf82d6472
function gini(x)
    n = length(x)
    xx = sort(x)
    2*(sum(collect(1:n).*xx))/(n*sum(xx))-1
end

function empirical_interaction_strength(B::Vector{Float64}, params::ModelParameters)

    S = size(params.network.species, 1)
    int = zeros(S, S)
    h = params.functional_response.h
    B0 = params.functional_response.B0
    ω = params.functional_response.ω

    for i in 1:S
        int[i, :] = [
                     ( params.biorates.x[i] * params.biorates.y[i] * B[i] * ω[i, j] * (B[j])^h ) /
                     ( (B0[i])^h + params.functional_response.c[i] * B[i] * (B0[i])^h  + sum(ω[i,:] .* (B .^h)))
                     for j in 1:S
                    ]
    end
    sparse(int)
end

function empirical_interaction_strength(solution, params::ModelParameters; last::Int64 = 1000)
    bm = foodweb_cv(ti, last = last).bm_sp

    empirical_interaction_strength(bm, params)

end
