using MultiScaleArrays
using DifferentialEquations
using LinearAlgebra
import Plots
import ForwardDiff
import Clustering
import LinearAlgebra
import RCall
import StatsBase


include("mu-types.jl")
include("mu-ecology.jl")
include("mu-evolution.jl")
include("general.jl")
R"library(bipartite)"
R"library(mvtnorm)"
R"source('R/toolbox_niche_difference_2022.R')"


plant1 = Plant([3.0])
plant2 = Plant([1.0])
plant3 = Plant([6.5])

animal1 = Animal([2.0])
animal2 = Animal([5.0])
animal3 = Animal([8.0])

pp = Population([plant1, plant2, plant3])
pa = Population([animal1, animal2, animal3])

tissue1 = Tissue(pp, pa) # Make a Tissue from Populations

pars =
    ParSetMu(
        sigma=1.0,
        a=1.0,
        r_plants= 1.0,
        r_animals= 0.0, 
        a_plants= -1.0,
        a_animals= -1.0,
        phi_plants= 0.01,
        phi_animals= 0.1,
        speed_a = 1.,
        speed_p = 1.)


get_r(tissue1, pars)
get_alpha(tissue1, pars)
nstar(tissue1, pars)

tissue1
w_i(animal1, tissue1, pars)
dw_i(Animal([3.0]), tissue1, pars)
dw_i(Plant([3.0]), tissue1, pars)


####

pp = Population([plant1, plant2, plant3])
pa = Population([animal1, animal2, animal3])
tissue1 = Tissue(pp, pa) # Make a Tissue from Populations

prob_mu = ODEProblem(evol1, tissue1, [0, 70.], Sp(deepcopy(pars)))
sol = solve(prob_mu, Tsit5(), callback = CallbackSet(singularstrategy(),extcallback), abstol=1e-14, reltol=1e-12, saveat = 0:1:5000, maxiters=1e7)

# sol_mu = solve(prob_mu, Tsit5(), callback=CallbackSet(extcallback, singularstrategy(), mubranch_callback), abstol=1e-14, reltol=1e-12, maxiters=1e7)
# [first.(get_plants(x).nodes) for x in sol.u]
# [first.(get_animals(x).nodes) for x in sol.u]

Plots.scatter(get_tidy(sol.u, sol.t, 1)..., color="green", label="plants")
Plots.scatter!(get_tidy(sol.u, sol.t, 2)..., color="orange", label="pollinators")
p1=Plots.plot!(xlabel = "time", ylabel = "trait value",legend=:bottomleft)
Plots.savefig("mu-dynamics.png")

# net = get_gamma(sol.u[end], pars)
# nestedness(kmeans_bin(net))
# global_stability(get_alpha(sol.u[end],pars))

alphas = [get_alpha(t, sol.prob.p.p) for t in sol.u]
gammas = [get_gamma(t, sol.prob.p.p) for t in sol.u]

p2=
Plots.scatter(sol.t, [OmegaL1.(alphas) bp_wine.(gammas) bp_modularity.(gammas) StatsBase.mean.(gammas)], 
                layout=(4,1), legend=:right, ylabel = ["Omega" "wine" "Q" "gamma"], 
                labels=["Structural stability" "W Nestedness" "Modularity" "Mutualistic strength"])

Plots.plot(p1, p2, layout=(1, 2), size=(1000,700))
Plots.savefig("mu-dynamics.png")

####

# # nstar(Population([Plant([10*rand()]) for _ in 1:10]), pars)
# remove_node!(tissue1, 2, 3)
# add_node!(tissue1,animal1,2)
# # remove_node!(tissue1, 2, 2)
# get_cart(tissue1, 6)
# remove_node!(tissue1,(1,2)...)
# # vcat([[(j,i) for i in 1:x] for (j,x) in pairs((3,9))]...)
# # [println(i,v) for (i,v) in pairs((2,8))]