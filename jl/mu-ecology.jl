######################################################################

# Ecology


"Area under the product of two gaussians"
function gaussians_overlap(m1::Real, m2::Real, s1::Real, s2::Real, a1::Real, a2::Real)
    a1 * a2 * sqrt(2) * s1 * s2 * exp(-(m1 - m2)^2 / (2 * s1^2 + 2 * s2^2)) * (1 / sqrt(s1^2 + s2^2))
end

alpha(x::Plant, y::Plant, p)::Real = p.a_plants * gaussians_overlap(x[1], y[1], p.sigma, p.sigma, p.a, p.a)
alpha(x::Animal, y::Animal, p)::Real = p.a_animals * gaussians_overlap(x[1], y[1], p.sigma, p.sigma, p.a, p.a)
alpha(x::Plant, y::Animal, p)::Real = p.phi_plants * gaussians_overlap(x[1], y[1], p.sigma, p.sigma, p.a, p.a)
alpha(x::Animal, y::Plant, p)::Real = p.phi_animals * gaussians_overlap(x[1], y[1], p.sigma, p.sigma, p.a, p.a)


rf(x::Plant, p) = 3 * gaussians_overlap(x[1], 0., p.sigma, 10., p.a, 1.)
# rf(x::Plant, p) = p.r_plants
rf(x::Animal, p) = p.r_animals

get_r(tissue::Tissue, p) = [rf(z1, p) for z1 in level_iter(tissue, 2)]
get_alpha(tissue::Tissue, p) = [alpha(z1, z2, p) for z1 in level_iter(tissue, 2), z2 in level_iter(tissue, 2)]

# Watch out for sign convention of alpha matrix !

nstar(t::Tissue, p)::Vector{Float64} = (-1 * get_alpha(t, p)) \ get_r(t, p)

# import Distributions
# "gaussian r"
# function r_normal(x, mu, sigma, rmax)
#     d = Distributions.Normal(mu, sigma)
#     (rmax / (Distributions.pdf(d, mu))) * Distributions.pdf(d, x)
# end