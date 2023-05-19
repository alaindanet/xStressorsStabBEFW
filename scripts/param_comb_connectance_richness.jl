using CSV, StatsBase, DataFrames, Arrow

rep = 1:50
S = [5, 10, 20, 40, 60]
C = 0.02:.05:.32
sigma = .2:.2:1.0
Z = [1, 2, 5, 10, 20, 40, 100]
ρ = 0:.2:1
K = [5, 10, 20, 30]

names = (:rep, :richness, :connectance, :Z, :sigma, :rho, :K)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(rep, S, C, Z, sigma, ρ, K)
           )[:]


# Filter impossible combination of C/S
limitCS = (
           S = S,
           Cmin = round.([(i - 1)/ i^2 + .01 for i in S], digits = 2),
           Cmax = [.32, .31, .24, .15, .09, .07]
          )
# Select good combinations
goodCSparam_v = [
                 findall(x ->
                         (
                          x.connectance >= limitCS.Cmin[i] &&
                          x.connectance <= limitCS.Cmax[i]) &&
                         x.richness == limitCS.S[i],
                         param
                        )
                 for i in 1:length(limitCS.S)
                ]
goodCSparam_idxs = reduce(vcat, goodCSparam_v)
bad_param = param[1:length(param) .∉ Ref(goodCSparam_idxs)]
good_param = param[
                   StatsBase.sample(goodCSparam_idxs,
                                    length(goodCSparam_idxs),
                                    replace=false
                                   )
                  ]

Arrow.write("../param_comb_connectance_richness.arrow", DataFrame(good_param))
ti = DataFrame(Arrow.Table("../param_comb_connectance_richness.arrow"))
