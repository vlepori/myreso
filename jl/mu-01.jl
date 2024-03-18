using MultiScaleArrays
using DifferentialEquations
import Plots
import ForwardDiff

include("mu-types.jl")
include("mu-ecology.jl")
include("mu-evolution.jl")

plant1 = Plant([3.0])
plant2 = Plant([1.0])
plant3 = Plant([6.5])

animal1 = Animal([2.0])
animal2 = Animal([5.0])
animal3 = Animal([8.0])

pp = Population([plant1, plant2, plant3])
pa = Population([animal1, animal2, animal3])

tissue1 = Tissue(pp, pa) # Make a Tissue from Populations

pars = Dict{Symbol,Any}()
push!(pars, :sigma => 1.0)
push!(pars, :a => 1.0)
push!(pars, :r_plants => 3.0)
push!(pars, :r_animals => 0.)
push!(pars, :ext => false)
push!(pars, :extmorph => 0)
push!(pars, :speed_a => 1.0)
push!(pars, :speed_p => 1.0)


get_r(tissue1, pars)
get_alpha(tissue1, pars)
nstar(tissue1, pars)

tissue1

w_i(animal1, tissue1, pars)
dw_i(Animal([3.0]), tissue1, pars)
dw_i(Plant([3.0]), tissue1, pars)
nstar(tissue1, pars)


####

pp = Population([plant1, plant2, plant3])
pa = Population([animal1, animal2, animal3])
tissue1 = Tissue(pp, pa) # Make a Tissue from Populations

prob_mu = ODEProblem(evol1, tissue1, [0, 40.], deepcopy(pars))
sol = solve(prob_mu, Tsit5(), callback = extcallback, abstol=1e-14, reltol=1e-12, saveat = 0:1:50, maxiters=1e7)

# sol_mu = solve(prob_mu, Tsit5(), callback=CallbackSet(extcallback, singularstrategy(), mubranch_callback), abstol=1e-14, reltol=1e-12, maxiters=1e7)
# [first.(get_plants(x).nodes) for x in sol.u]
# [first.(get_animals(x).nodes) for x in sol.u]

Plots.scatter(get_tidy(sol.u, sol.t, 1)..., color="green", label="plants")
Plots.scatter!(get_tidy(sol.u, sol.t, 2)..., color="orange", label="pollinators")
Plots.plot!(xlabel = "time", ylabel = "trait value",legend=:bottomleft)
Plots.savefig("mu-dynamics.png")


####

# nstar(Population([Plant([10*rand()]) for _ in 1:10]), pars)
remove_node!(tissue1, 2, 3)
add_node!(tissue1,animal1,2)
# remove_node!(tissue1, 2, 2)
get_cart(tissue1, 6)
remove_node!(tissue1,(1,2)...)
# vcat([[(j,i) for i in 1:x] for (j,x) in pairs((3,9))]...)
# [println(i,v) for (i,v) in pairs((2,8))]