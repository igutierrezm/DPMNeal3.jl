module DPMNeal3

using Parameters: @with_kw, @unpack
using Random: AbstractRNG, randperm!

# @with_kw struct DPMModel
#     N::Int
#     K::Ref{Int} = Ref(1)
#     α::Ref{Float64} = Ref(1.0)
#     τ::Vector{Int} = collect(1:N)
#     d::Vector{Int} = ones(N)
#     n::Vector{Int} = [N]
#     D::Set{Int} = Set(1)
# end

# # Accessors

# N(m::DPMModel) = m.N
# d(m::DPMModel) = m.d
# n(m::DPMModel) = m.n
# D(m::DPMModel) = m.D
# α(m::DPMModel) = m.α

# # Interface

# logpdf0(m::DPMModel, y, i) = error("Not implemented")
# logpdf1(m::DPMModel, y, i, j) = error("Not implemented")

# function update_τ!(m::DPMModel, rng)
#     randperm!(rng, m.τ)
# end

# function update_d!(m::DPMModel, rng, y)
#     @unpack D, d, n, α = m
#     for i = 1:N0
#         d0 = d[i]
#         p0 = -Inf
#         k0 = 0
#         for k in D
#             p1 = logpdf1(m, y, i, j) + log(n[k] - (d0 == k))
#             p1 = p1 - log(-log(rand(rng)))
#             p1 < p0 && continue
#             d0[i] = k
#             p0 = p1
#         end
#         p1 = logpdf1(m, y, i, k0) + log(α[])
#         p1 = p1 - log(-log(rand(rng)))
#         p1 > p0 && (d[i] = k0)
#     end
# end

# function update!(m::DPMModel, rng, y)
#     update_τ!(m, rng)
#     update_d!(m, rng, y)
# end

end # module
