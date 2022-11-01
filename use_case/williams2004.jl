using BEFWM2
using DifferentialEquations
using Plots
include("src/minmax.jl")

Random.seed!(1234)

function connectance(mat)
    sum(mat) / size(mat,1)^2
end

foodweb = FoodWeb(nichemodel, 10,C = .14)
while connectance(foodweb.A) != .14 
    global foodweb = FoodWeb(nichemodel, 10,C = .14)
end
# Generally, reltol is the relative accuracy while abstol is the accuracy when u
# is near zero. These tolerances are local tolerances and thus are not global
# guarantees. However, a good rule of thumb is that the total solution
# accuracy is 1-2 digits less than the relative tolerances. Thus for the
# defaults abstol=1e-6 and reltol=1e-3, you can expect a global accuracy
# of about 1-2 digits.


params = ModelParameters(foodweb, functional_response = BioenergeticResponse(foodweb, h = 0.025 + 1))
B = rand(size(foodweb.A, 1))
out = []
q0 = 0.0:.05:3
for q in q0
    println("Start q = $q, h = $(q+1)")
    params_tmp = ModelParameters(foodweb, functional_response = BioenergeticResponse(foodweb, h = q + 1))
    sol = simulate(
                   params_tmp,
                   rand(size(foodweb.A, 1)),# repeat([.1], size(foodweb.A, 1)),
                   alg = Vern9(),
                   tmax = 2000,
                   callback = CallbackSet(
                                          BEFWM2.PositiveDomain(),
                                          # TerminateSteadyState(1e-6, 1e-4),
                                          BEFWM2.ExtinguishSpecies(1e-10, false),
                                         ),
                   # alg_hints = [:stiff], 
                   # reltol=1e-8, abstol=1e-10
                  )
    tmp = []
    for i in 1:size(foodweb.A, 1)
        push!(tmp, (sp = i, extrema = unique(local_extrema(sol, sp = i, nt = 10000, tmin = Int(round(maximum(sol.t) - 1000)), tmax = Int(round(maximum(sol.t))), verbose = false).extrema)));
    end
    push!(out, (q0 = q, sp_extrema = tmp))
    println("End q = $q, h = $(q+1)")
end

p = []
for s in 1:size(foodweb.A, 1)
    local tmp = scatter(
                  repeat([out[1].q0], length(out[1].sp_extrema[s].extrema)),
                  out[1].sp_extrema[s].extrema,
                  legend = false,
                  xlabel = "q",
                  ylabel = "Bmin, Bmax"
                 )
    push!(p, tmp)
    for i in 2:length(out)
        p[s] = scatter!(
                      repeat([out[i].q0], length(out[i].sp_extrema[s].extrema)),
                      out[i].sp_extrema[s].extrema,
                      legend = false)
    end
end
length(p)
plot(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])

