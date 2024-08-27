# ######################################################################


Base.@kwdef struct ParSetFullB <: AbstractParSet
    # Competition
    sigma_comp_plants::Tuple{Float64,Float64}
    sigma_comp_animals::Tuple{Float64,Float64}
    a_comp_plants::Float64
    a_comp_animals::Float64

    # Mutualism
    sigma_mutualism::Float64
    a_mutualism_plants::Float64
    a_mutualism_animals::Float64
    d::Float64

    # Rs
    mu_r_plants::Tuple{Float64,Float64}
    mu_r_animals::Tuple{Float64,Float64}
    sigma_r_plants::Tuple{Float64,Float64}
    sigma_r_animals::Tuple{Float64,Float64}
    r_max_plants::Float64
    r_max_animals::Float64

    # Balance
    delta::Float64
    speed_a::Float64
    speed_p::Float64
end




# Ecology under model 2: symmetric - asymmetric interactions

"""
    difference(x, y, d)

Values close to 0 for x >> y, 1 for x << y. Sharp for d -> 0.
"""
function difference(x, y, d)
    1/(1+exp(d*(x - y)))
end

"Gaussian kernel, max = 1"
function matching(m1::Real, m2::Real, a::Real, s::Real)
   a * exp(-(m1 - m2)^2 / (2s^2))
end

function alpha(x::Plant, y::Plant, p)::Real
    ax1 = matching(x[1], y[1], 1, p.sigma_comp_plants[1])
    ax2 = matching(x[2], y[2], 1, p.sigma_comp_plants[2])

    p.a_comp_plants * wgmean(ax1, ax2, 1, p.delta)
end

function alpha(x::Animal, y::Animal, p)::Real
    ax1 = matching(x[1], y[1], 1, p.sigma_comp_animals[1])
    ax2 = matching(x[2], y[2], 1, p.sigma_comp_animals[2])

    p.a_comp_animals * wgmean(ax1, ax2, 1, p.delta)
end


function alpha(x::Animal, y::Plant, p)::Real
    ax1 = matching(x[1], y[1], 1, p.sigma_mutualism)
    ax2 = difference(y[2], x[2], p.d)

    p.a_mutualism_animals * wgmean(ax1, ax2, 1, p.delta)
end

function alpha(x::Plant, y::Animal, p)::Real
    ax1 = matching(x[1], y[1], 1, p.sigma_mutualism)
    ax2 = difference(x[2], y[2], p.d)

    p.a_mutualism_plants * wgmean(ax1, ax2, 1, p.delta)
end


function rf(x::Plant, p)::Real
    ax1 = matching(x[1], p.mu_r_plants[1], 1, p.sigma_r_plants[1])
    ax2 = matching(x[2], p.mu_r_plants[2], 1, p.sigma_r_plants[2])

    p.r_max_plants * wgmean(ax1, ax2, 1, p.delta)
end


function rf(x::Animal, p)::Real
    ax1 = matching(x[1], p.mu_r_animals[1], 1, p.sigma_r_animals[1])
    ax2 = matching(x[2], p.mu_r_animals[2], 1, p.sigma_r_animals[2])

    p.r_max_animals * wgmean(ax1, ax2, 1, p.delta)
end






# function alpha(x::Animal, y::Animal, p)::Real
#     ax1 = matching(x[1], y[1], p.a_competition_animals, p.sigma_competition)
#     ax2 = difference(x[2], y[2], p.d)
#     wg_mean(ax1, ax2, 1, p.delta)
# end

# function alpha(x::Plant, y::Animal, p)::Real
#     ax1 = matching(x[1], y[1], p.a_competition_plants, p.sigma_competition)
#     ax2 = difference(x[2], y[2], p.d)
#     wg_mean(ax1, ax2, 1, p.delta)
# end

# function alpha(x::Animal, y::Plant, p)::Real
#     ax1 = matching(x[1], y[1], p.a_competition_plants, p.sigma_competition)
#     ax2 = difference(x[2], y[2], p.d)
#     wg_mean(ax1, ax2, 1, p.delta)
# end




# alpha(x::Animal, y::Animal, p)::Real = matching(x[1], y[1], p.a_competition_animals, p.sigma_competition)
# alpha(x::Plant, y::Animal, p)::Real = matching(x[1], y[1], p.a_mutualism_plants, p.sigma_mutualism)
# alpha(x::Animal, y::Plant, p)::Real = matching(x[1], y[1], p.a_mutualism_animals, p.sigma_mutualism)


# rf(x::Plant, p) = 3 * gaussians_overlap(x[1], 0.0, p.sigma_plants, 100.0, p.a, 1.0)




# wgmean(5, 0, 10, 1)
# wgmean(x::Vector) = wgeomean(x, ones(length(x)))
# wgmean([5, 1], [10, 100])





# # B(x) = exp( −1/(1−x²) )

# function alpha(x::Plant, y::Plant, p)::Real
#     -sym_kernel(x[1],y[1], 2.)
# end

# function alpha(x::Animal, y::Animal, p)::Real
#     -sym_kernel(x[1], y[1], 1.0)*0.1
# end

# function alpha(x::Plant, y::Animal, p)::Real
#     sym_kernel(x[1], y[1], 1.0)*asym_kernel(x[2],y[2],5.)
# end

# function alpha(x::Animal, y::Plant, p)::Real
#     sym_kernel(x[1], y[1], 1.0)*asym_kernel(x[2], y[2], 5.0)
# end


# function rf(x::Plant, p)
#     5*sym_kernel(x[1],0.,10)*sym_kernel(x[2],0.,10)
# end
# # rf(x::Plant, p) = p.r_plants
# rf(x::Animal, p) = p.r_animals

# # get_r(tissue::Tissue, p) = [rf(z1, p) for z1 in level_iter(tissue, 2)]
# # get_alpha(tissue::Tissue, p) = [alpha(z1, z2, p) for z1 in level_iter(tissue, 2), z2 in level_iter(tissue, 2)]

# # Watch out for sign convention of alpha matrix !

# # nstar(t::Tissue, p)::Vector{Float64} = (-1 * get_alpha(t, p)) \ get_r(t, p)
# # nstar(alpha::Matrix, r::Vector) = (-1 * alpha) \ r

# # global_stability(alpha::Array) = minimum(LinearAlgebra.eigen(-1 * alpha).values) > 0 ? 1 : 0
# # global_stability(t::Tissue, p) = global_stability(get_alpha(t,p))

# # function lv(du, u, p, t)
# #     r,alpha = p
# #     du[:] = u .* (r + alpha * u)
# # end

# # function ecol_dynamics(tissue, p)
# #     tspan = (0.0, 50.)
# #     paramet = (get_r(tissue, p), get_alpha(tissue, p))
# #     prob = ODEProblem(lv, 5*rand(length(paramet[1])), tspan, paramet)
# #     solve(prob)
# # end

# # Plots.plot(ecol_dynamics(tissue,pars))
# # nstar(tissue1, pars)

# # # import Distributions
# # # "gaussian r"
# # # function r_normal(x, mu, sigma, rmax)
# # #     d = Distributions.Normal(mu, sigma)
# # #     (rmax / (Distributions.pdf(d, mu))) * Distributions.pdf(d, x)
# # # end