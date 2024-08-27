# using MultiScaleArrays

function check_values(v::Vector{T} where T <: Real)::Bool
    (length(v) == 2) # & (last(v) > 0.)
end

struct Plant{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    Plant(values::Vector{S}) where {S<:Real} = check_values(values) ? new{S}(values) : error("Define two values")
end

struct Animal{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    Animal(values::Vector{S}) where {S<:Real} = check_values(values) ? new{S}(values) : error("Define two values")
end

struct Population{T <: AbstractMultiScaleArrayLeaf, B <: Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

Population(v::Vector{T} where {T<:AbstractMultiScaleArrayLeaf}) = length(unique(typeof.(v))) == 1 ? construct(Population, deepcopy(v)) : error("1 type per population")


struct Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

# Safer constructor

function Tissue(p1::Population, p2::Population)
    if prod(is_plant.(get_nodes(p1))) & prod(is_animal.(get_nodes(p2)))
        construct(Tissue, deepcopy([p1, p2])) # Make a Tissue from Populations
    else
        error("Provide one plant and one animal population")
    end
end

function Tissue(vp::Vector{Vector{Float64}}, va::Vector{Vector{Float64}})
    Tissue(
        Population([Plant(x) for x in vp]),
        Population([Animal(x) for x in va])
    )
end

get_plants(t::Tissue) = t.nodes[1]
get_animals(t::Tissue) = t.nodes[2]

get_nodes(p::Population) = p.nodes
is_animal(x) = isa(x, Animal)
is_plant(x) = isa(x, Plant)



# Pretty printing:

import Base.show

function Base.show(io::IO, c::Population)
    println(io, "Population, size = $(length(c.nodes)).")
    for n in c.nodes
        println(io, "    $n")
    end
end

function Base.show(io::IO, c::Tissue)
    println(io, "Tissue of $(length(c.nodes)) populations.\n")
    for n in c.nodes
        println(io, "    $n")
    end
end

# function Base.show(io::IO, v::Vector{T} where {T<:Tissue})
function Base.show(io::IO, v::Vector{Tissue{Population{T,Float64} where T<:AbstractMultiScaleArrayLeaf,Float64}})
    for n in v
        println(io, "Tissue, values: $(n[:])")
    end
end

abstract type AbstractParSet end

# Base.@kwdef struct ParSetFull <: AbstractParSet
#     # sigma::Float64
#     # a::Float64
#     delta::Float64
#     r_plants::Float64
#     r_animals::Float64
#     a_plants::Float64
#     a_animals::Float64
#     phi_plants::Float64
#     phi_animals::Float64
#     m_plants::Float64
#     m_animals::Float64
#     speed_a::Float64
#     speed_p::Float64
# end












# ######################################################################

# plant1 = Plant([5.0, 4.0])
# plant2 = Plant([2.0, 1.0])
# # plant1 = Plant([5.0, 3.0])

# animal1 = Animal([5.0, 3.0])
# animal2 = Animal([2.0, 1.0])
# # plant2 = Plant([2.0, 1.0])

# pp = construct(Population, deepcopy([plant1, plant2]))
# pa = construct(Population, deepcopy([animal1, animal2]))

# # tissue1
# tissue1 = construct(Tissue, deepcopy([pp, pa])) # Make a Tissue from Populations

# f(x::Plant) = println("Hey I'm a plant")
# f(x::Animal) = println("BAUUU")

# rf(x::Plant) = 5
# rf(x::Animal) = 2


# alpha(x::Plant, y::Plant)::Real = 1
# alpha(x::Animal, y::Animal)::Real = 2
# alpha(x::Plant, y::Animal)::Real = 0
# alpha(x::Animal, y::Plant)::Real = 0.5

# get_alpha(tissue::Tissue) = [alpha(z1, z2) for z1 in level_iter(tissue, 2), z2 in level_iter(tissue, 2)]


# Tissue(pa,pp)

# # Iterate at the second level of tissue to find the organisms
# for x in level_iter(tissue1, 2)
#     println(typeof(x))
#     f(x)
#     # Do something with the cells!
# end

# [rf(z) for z in level_iter(tissue1, 2)]
# # [alpha(z1, z2) for z1 in level_iter(tissue1.nodes[1], 1), z2 in level_iter(tissue1.nodes[2], 1)]


# # function nstar(x::Tissue, p)

# #     # n = length(x)
# #     np = p[:np] # of plants
# #     w, d, h = p[:w], p[:d], p[:h]
# #     asig, amax = p[:asig], p[:amax]
# #     rmu, rsig, rmax = p[:rmu], p[:rsig], p[:rmax]
# #     xp = x[1:np]
# #     xa = x[(np+1):end]
# #     rp = r_normal.(xp, rmu, rsig, rmax)
# #     ra = fill(0.1, length(x) - np)
# #     # ra = p[:ra]
# #     r = [rp; ra]

# #     network_strength = [olap(xi, xj, w, d, h) for xi in xp, xj in xa] #|> x-> x/maximum(x)
# #     ## benefits to plants
# #     gamma_p = 0.05 * network_strength
# #     ## benefits to animals
# #     gamma_a = 0.01 * copy(network_strength')

# #     alpha_p = [alpha_normal(xi, xj, asig, amax) for xi in xp, xj in xp]
# #     # alpha_a = p[:alpha_a]
# #     alpha_a = [olap(xi, xj, w, d, h) for xi in xa, xj in xa]
# #     alpha = make_alpha(alpha_p, alpha_a, gamma_p, gamma_a)

# #     # return eigen(alpha).values
# #     min_eig = minimum(real.(eigen(alpha).values))
# #     # println(min_eig)
# #     # min_eig>0|| error("instability")
# #     N = alpha \ r
# #     # minimum(N) > 0 || @warn (minimum(N),"unfeasibility!")

# #     return N
# # end
# # cell3 = Morph([3.0; 2.0; 5.0])
# # cell4 = Morph([4.0; 6.0])
# # population = construct(Population, deepcopy([cell1, cell3, cell4]))
# # population2 = construct(Population, deepcopy([cell1, cell3, cell4]))
# # population3 = construct(Population, deepcopy([cell1, cell3, cell4]))
# # tissue1 = construct(Tissue, deepcopy([population, population2, population3])) # Make a Tissue from Populations
# # tissue2 = construct(Tissue, deepcopy([population2, population, population3]))
