########################################

# Evolution

function w_i(m, t, N::Vector, p)
    rf(m, p) + sum([alpha(m, z1, p) for z1 in level_iter(t, 2)] .* N)
end

w_i(m, t::Tissue, p) = w_i(m, t, nstar(t, p), p)

# ??? define methods twice for plants and animals?
dw_i(m::Animal, t, N, p) = ForwardDiff.gradient(x -> w_i(Animal(x), t, N, p), [m[1]])
dw_i(m::Animal, t::Tissue, p) = dw_i(m, t, nstar(t, p), p)

dw_i(m::Plant, t, N, p) = ForwardDiff.gradient(x -> w_i(Plant(x), t, N, p), [m[1]])
dw_i(m::Plant, t::Tissue, p) = dw_i(m, t, nstar(t, p), p)


function evol1(du, u, p, t)
    N = nstar(u, p)

    if any(N .< 0.0001)
        p[:ext] = true
        extidx = findall(x -> x < 0.0001, N)
        length(extidx) > 1 && error("multiple extinctions!")
        p[:extmorph] = first(extidx)
        # error("Extinction at $(p[:extmorph]).")
    end

    for (m, dm) in LevelIter(2, u, du)
        dm.values[:] = dw_i(m, u, N, p)
    end

end

# Plotting trajectories

"level in [1,2]. returns two vectors"
function get_tidy(u, t, level::Int)
    out = [[0.0 0.0]]
    for i in 1:length(t)
        for n in u[i].nodes[level].nodes
            push!(out, [t[i] n[1]])
        end
    end
    out2 = vcat(out[2:end, :]...)
    return out2[:, 1], out2[:, 2]
end

# Get indexes for deleting from linear indices. t is a tissue or integrator.u
"Return a 2-Tuple"
get_cart(t)::Vector{Tuple{Int64,Int64}} = vcat([[(j, i) for i in 1:x] for (j, x) in pairs(t.end_idxs)]...)

get_cart(t, i::Int)::Tuple{Int64,Int64} = get_cart(t)[i]



# Extinction callbacks

extcondition(u, t, integrator) = integrator.p[:ext] # or make continuous?

function extaffect!(integrator)
    
    rmidx = get_cart(integrator.u, integrator.p[:extmorph])
    println("EXTINCTION: $rmidx")
    remove_node!(integrator, rmidx...)
    
    integrator.p[:extmorph] = 0
    integrator.p[:ext] = false
    # printnice(integrator.u)
end

extcallback = DiscreteCallback(extcondition, extaffect!)


# d2w_i(m, t, N, p) = ForwardDiff.hessian(x -> w_i(Morph(x[1], m[2]), t, N, p), [m[1]])
# dw_i(m, t, N, p) = ForwardDiff.gradient(x -> w_i(Morph(x[1], m[2]), t, N, p), [m[1]])
# d2w_i(m, t, N, p) = ForwardDiff.hessian(x -> w_i(Morph(x[1], m[2]), t, N, p), [m[1]])
