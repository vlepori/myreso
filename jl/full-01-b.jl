# This is the hybrid model (symmetric-asymmetric)

using MultiScaleArrays
using DifferentialEquations
using RCall
using DataFrames


import LinearAlgebra
import Plots
import ForwardDiff
import Clustering
import Random
import LinearAlgebra
import StatsBase

include("full-types.jl")
include("full-ecology-b.jl")
include("full-evolution.jl")
include("general.jl")
# R"library(bipartite)"
# R"library(mvtnorm)"
# R"library(pracma)"
# R"source('R/toolbox_coexistence_2024.R')"

plant1 = Plant([-3.0,1.5])
plant2 = Plant([2.0,0.5])
plant3 = Plant([0.5, 2.5])

animal1 = Animal([0.0,5.])
animal2 = Animal([5.0, 2.])
animal3 = Animal([3.0, 1.5])

pp = Population([plant1,  plant3])
pa = Population([animal1, animal3])

tissue1 = Tissue(pp, pa) # Make a Tissue from Populations


pars = ParSetFullB(
        # Competition
        sigma_comp_plants = (0.8, 0.8),
        sigma_comp_animals = (0.8, 0.8),
        a_comp_plants = -1.3,
        a_comp_animals = -.4,
        # Mutualism
        sigma_mutualism = 1.0,
        a_mutualism_plants  = 0.2, # mutualisic gain plants
        a_mutualism_animals = 0.3,
        d = 1., # steepness of asymmetric function
        # Rs
        mu_r_plants = (0.0, 3.0),
        mu_r_animals = (0.0, -3.0),
        sigma_r_plants = (1.0, 1.0),
        sigma_r_animals = (1.0, 1.0),
        r_max_plants = 5.0,
        r_max_animals = 1.0,
        # Balnce
        delta = 1, # asymmetric strength, sym strenght set to 1
        speed_a = 1.,
        speed_p = 1.
)

get_alpha(tissue1, pars)
get_gamma(tissue1, pars)
get_r(tissue1, pars)
nstar(tissue1, pars)
rmorph(rangemu::Array, rangesigma::Array) = [runif(rangemu), runif(rangesigma)]
rmorph(n::Int, rangemu, rangesigma) = [rmorph(rangemu, rangesigma) for _ in 1:n]
Plots.heatmap(get_gamma(Tissue(rmorph(5, [-2, 2],[0,5]), rmorph(5, [-2, 2],[0,.5])), pars))
w_i(animal1, tissue1, pars)
dw_i(Plant([3.0, 0.5]), tissue1, pars)

######################################################

plant1 = Plant([1.0, 0.2])
plant2 = Plant([-1.0, 0.8])
animal1 = Animal([-0.5, 0.3])
animal2 = Animal([-0.8, 5.5])
pp = Population([plant1, plant2])
pa = Population([animal1,animal2])
tissue = Tissue(pp, pa) # Make a Tissue from Populations
nstar(tissue,pars)

prob_mu = ODEProblem(evol, tissue, [0, 50.], Sp(deepcopy(pars)))
sol = solve(prob_mu, Tsit5(), saveat = 0:5:2500, abstol = 1e-14, reltol=1e-12, maxiters=1e7)

plot_p = get_tidy(sol.u, sol.t, 1)
plot_a = get_tidy(sol.u, sol.t, 2)
Plots.scatter(plot_p[2], plot_p[3], markerstrokewidth=1,color="green", label="plants")
Plots.scatter!(plot_a[2], plot_a[3], markerstrokewidth=.5, color="orange", label="pollinators")

# Plots.heatmap(get_gamma(sol.u[1], pars))
# Plots.heatmap(get_gamma(sol.u[end], pars))

#########################################################

xt = Tissue(sort(rmorph(50, [-5, 5], [-5, 5])), sort(rmorph(50, [-5, 5], [-5, 5])))

Plots.heatmap(get_gamma(sort_x1(xt), pars))
Plots.heatmap(get_gamma(sort_x2(xt), pars))


##########################################################

autoheatmap(-5:0.1:5, -5:0.1:5, (x, y) -> rf(Plant([x, y]), pars))
autoheatmap(-5:0.1:4, -5:0.1:5, (x, y) -> rf(Animal([x, y]), pars))

sum([is_feasible(Tissue(rmorph(7,[-5,5], [-5, 5]), rmorph(7, [-5,5],[-5, 5])), pars) for _ in 1:100])

folder = "sims-full-mix/30aug24/"
mkdir(folder)


for r in 1:10000

        np = 4
        na = 4

        tissue0 = Tissue(sort(rmorph(np, [-5, 5], [-5, 5])), sort(rmorph(na, [-5, 5], [-5, 5])))
        id = folder * lpad(r, 4, '0') * "-" * Random.randstring()

        if is_feasible(tissue0, pars)

                try
                        pr = ODEProblem(evol, tissue0, [0, 100.0], Sp(deepcopy(pars)))
                        sl = solve(pr, Tsit5(), callback=extcallback, abstol=1e-12, reltol=1e-10, saveat=0:1:5000, maxiters=1e7)

                        ex = richness(sl.u[1]) - richness(sl.u[end])

                        save_sol(sl.u[1], sl.prob.p.p, sl.t[1], ex, id * "-s.rds")
                        save_sol(sl.u[end], sl.prob.p.p, sl.t[end], ex, id * "-f.rds")

                        println("saved $id, r = $r")

                catch e
                        println(e)
                end

        end
end
