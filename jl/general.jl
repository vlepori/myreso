# ECOLOGY; ANY MODEL 

get_r(tissue::Tissue, p) = [rf(z1, p) for z1 in level_iter(tissue, 2)]
get_alpha(tissue::Tissue, p) = [alpha(z1, z2, p) for z1 in level_iter(tissue, 2), z2 in level_iter(tissue, 2)]

get_gamma(t::Tissue, p) = [alpha(z1, z2, p) for z1 in level_iter(t.nodes[1], 1), z2 in level_iter(t.nodes[2], 1)]
# mean_gamma(tissue, p) = StatsBase.mean(get_gamma(tissue, p))

# Watch out for sign convention of alpha matrix !

nstar(t::Tissue, p)::Vector{Float64} = (-1 * get_alpha(t, p)) \ get_r(t, p)
nstar(alpha::Matrix, r::Vector) = (-1 * alpha) \ r

global_stability(alpha::Array) = minimum(LinearAlgebra.eigen(-1 * alpha).values) > 0 ? 1 : 0
global_stability(t::Tissue, p) = global_stability(get_alpha(t, p))

function lv(du, u, p, t)
    r, alpha = p
    du[:] = u .* (r + alpha * u)
end

function ecol_dynamics(tissue, p)
    tspan = (0.0, 50.0)
    paramet = (get_r(tissue, p), get_alpha(tissue, p))
    prob = ODEProblem(lv, 5 * rand(length(paramet[1])), tspan, paramet)
    solve(prob)
end

# NETWORK METRICS ######################################

"make a network of Int ∈ {0,1} from a continuous network"
function binarize(network)
    Int.(network .!= 0)
end

"count shared interactions in binary arrays"
function match_hits(A::Array{Int}, B::Array{Int})
    sum(A[A.==B])
end

"make two matrices (one per dimension) of shared interaction counts from binary network"
function shared_interactions(network::Array{Int})
    np, na = size(network)
    m_p = [match_hits(network[i, :], network[j, :]) for i in 1:np, j in 1:np]
    m_a = [match_hits(network[:, i], network[:, j]) for i in 1:na, j in 1:na]
    return m_p, m_a
end

"Sum the elements on the upper triangle of a square matrix `M`"
function sum_tri(M::Array)
    m, n = size(M)
    m == n || error("m not square")
    out = 0
    for i in 1:(n-1)
        out += sum(M[i, (i+1):end])
    end
    return out
end

"make pairwise matrices of min(links_i, links_j)"
function min_links(network::Array{Int})
    np, na = size(network)
    sum_p = sum(network, dims=2)
    sum_a = sum(network, dims=1)
    out_p = [min(sum_p[i], sum_p[j]) for i in 1:np, j in 1:np]
    out_a = [min(sum_a[i], sum_a[j]) for i in 1:na, j in 1:na]
    return out_p, out_a
end

"Calculate the nestedness ∈ [0,1] for a discrete network, following Bastolla et al. 2009"
function nestedness(network::Matrix{Int})

    si_p, si_a = shared_interactions(network)
    ml_p, ml_a = min_links(network)

    num = sum_tri(si_p) + sum_tri(si_a)
    den = sum_tri(ml_p) + sum_tri(ml_a)

    return num / den
end

"Discretize (0,1) using kmeans. Make sure that there are both links and missing links !"
function kmeans_bin(m::Matrix)
    out = zeros(Int, size(m))
    clusters = Clustering.kmeans(vec(m)', 2)
    ctrs = clusters.centers .== maximum(clusters.centers)
    println("<k-means (k=2), centers are $(clusters.centers) >")
    out[:, :] = [ctrs[i] for i in clusters.assignments]
    return out
end

function bp_wnodf(m::Matrix{Float64})
    out = R"nested($m, method = 'weighted NODF')"
    convert(Float64,out)
end

function bp_wine(m::Matrix{Float64})
    out = R"wine($m)[['win']]"
    convert(Float64, out)
end

function bp_modularity(m::Matrix{Float64})
    out = R"computeModules($m)@likelihood"
    convert(Float64, out)
end

function omega_l1(alpha::Matrix{Float64})
    out=R"(Omega_L1($alpha))"
    convert(Float64, out)
end

# function eta(r::Vector{Float64}, alpha::Matrix{Float64})
#     out = R"(angle_to_border($r, $alpha))"
#     convert(Float64, out)
# end

function OmegaL1(alpha::Array)
    size(alpha)[1] == 1 && return missing
    alpha_n = alpha * diagm(vec(1 ./ sum(alpha, dims=1)))
    out = missing
    try
        out = log(det(alpha_n))
    catch
        print("!")
    end
    return out / log(10)
end


####################################################

# MISCELLANEOUS AND PLOTTING

"heatmap of f(x,y) over range xy"
function autoheatmap(range_x, range_y, f::Function)
    Plots.heatmap(range_x, range_y, [f(x, y) for x in range_x, y in range_y]')
end

function autocontour(range_x, range_y, f::Function, levels=[0])
    Plots.contour(range_x, range_y, [f(x, y) for x in range_x, y in range_y]', levels=levels, color=:black)
end

function autocontour!(range_x, range_y, f::Function, levels=[0])
    Plots.contour!(range_x, range_y, [f(x, y) for x in range_x, y in range_y]', levels=levels, color=:black)
end


"from plantypos.jl"
function vectorf!(xs, ys, dxs::Array{Float64,2}, dys::Array{Float64,2})
    x, y = meshgrid(xs, ys)
    z1 = vcat(dxs...)
    z2 = vcat(dys...)
    q1 = fill(0.0, length(z1))
    q2 = deepcopy(q1)
    for i in 1:length(z1)
        q1[i], q2[i] = [z1[i], z2[i]] |> x -> 0.002 * x / sqrt(sum(x .^ 2))
    end
    Plots.quiver!(x, y, quiver=(q1, q2), color=:black)
end


meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))


"Consistent interface with plotting functions in this project. Here, f should be a function from R2 -> R2"
function vectorf!(range_x, range_y, f::Function)
    m = [f(x, y) for x in range_x, y in range_y]
    dx = first.(m)
    dy = second.(m)
    vectorf!(range_x, range_y, dx, dy)
end

second(x) = x[2]




# gamma(x, y, p)::Union{Real,Missing} = missing
# gamma(x::Plant, y::Animal, p)::Real = alpha(x, y, p)
# get_gamma(tissue::Tissue, p) = [gamma(z1, z2, p) for z1 in level_iter(tissue, 2), z2 in level_iter(tissue, 2)]
