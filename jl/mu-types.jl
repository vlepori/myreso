struct Plant{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    Plant(values::Vector{S}) where {S<:Real} = length(values) == 1 ? new{S}(values) : error("Define one mu value")
end

struct Animal{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    Animal(values::Vector{S}) where {S<:Real} = length(values) == 1 ? new{S}(values) : error("Define one mu value")
end

struct Population{T<:AbstractMultiScaleArrayLeaf,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

Population(v::Vector{T} where {T<:AbstractMultiScaleArrayLeaf}) = length(unique(typeof.(v))) == 1 && construct(Population, deepcopy(v))


# Do we need Tissues at all ?
struct Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end


# Convenience

get_plants(t::Tissue) = t.nodes[1]
get_animals(t::Tissue) = t.nodes[2]

get_nodes(p::Population) = p.nodes
is_animal(x) = isa(x, Animal)
is_plant(x) = isa(x, Plant)


# Safer constructor

function Tissue(p1::Population, p2::Population)
    if prod(is_plant.(get_nodes(p1))) & prod(is_animal.(get_nodes(p2)))
        construct(Tissue, deepcopy([p1, p2])) # Make a Tissue from Populations
    else
        error("Provide one plant and one animal population")
    end
end


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

Base.@kwdef struct ParSetMu <: AbstractParSet
    sigma::Float64
    a::Float64
    r_plants::Float64
    r_animals::Float64
    # ext::Bool
    # extmorph::Int
    speed_a::Float64
    speed_p::Float64
end




# [first.(get_plants(x).nodes,2) for x in sol_mu.u]


