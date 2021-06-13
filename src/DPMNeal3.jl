module DPMNeal3

abstract type AbstractDPM end
abstract type AbstractDPM_SB end
const DPMSB = AbstractDPM_SB

using Distributions: Beta, Gamma 
using Parameters: @unpack
using Lazy: @forward
using Random: randperm, randperm!, AbstractRNG
using SpecialFunctions: loggamma
export AbstractDPM, DPM, NormalDPM, cluster_labels, cluster_sizes, n_clusters, dp_mass, active_clusters, passive_clusters, max_cluster_label, @forward_dpm_methods, update!

"""
    DPM(rng::AbstractRNG, N; K0 = 1, a0 = 2.0, b0 = 4.0)

Initialize a generic DPM with `N` observations, `K0` initial clusters and a 
`Gamma(a0, b0)` prior distribution for the mass parameter.
"""
struct DPM <: AbstractDPM
    N::Int             # sample size
    K::Vector{Int}     # number of clusters
    Q::Vector{Int}     # largest cluster label
    α::Vector{Float64} # DP mass parameter
    τ::Vector{Int}     # ordering of the observations
    d::Vector{Int}     # cluster labels
    n::Vector{Int}     # cluster sizes
    P::Set{Int}        # passive clusters
    A::Set{Int}        # active clusters
    a0::Float64        # shape parameter in α's prior
    b0::Float64        # rate parameter in α's prior
    function DPM(rng, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)
        @assert N >= K0
        @assert K0 >= 0
        @assert a0 >= 0
        @assert b0 >= 0
        K = [K0]
        Q = [K0]
        α = [1.0]
        τ = randperm(rng, N)
        d = rand(rng, 1:K0, N); d[1:K0] = 1:K0;
        n = [sum(d .== k) for k in 1:(K0 + 1)]
        P = Set(K0 + 1)
        A = Set(1:K0)
        new(N, K, Q, α, τ, d, n, P, A, a0, b0)
    end
end

"""
    cluster_labels(m::AbstractDPM)

Return the cluster labels.
"""
cluster_labels(m::DPM) = m.d

"""
    cluster_sizes(m::AbstractDPM)

Return the cluster sizes.
"""
cluster_sizes(m::DPM) = m.n

"""
    active_clusters(m::AbstractDPM)

Return the set of active clusters.
"""
active_clusters(m::DPM) = m.A

"""
    passive_clusters(m::AbstractDPM)

Return a subset of the passive clusters.
"""
passive_clusters(m::AbstractDPM) = m.P

"""
    n_clusters(m::DPM)

Return the number of active clusters.
"""
n_clusters(m::DPM) = m.K[1]

"""
    max_cluster_label(m::AbstractDPM)

Return the largest (active) cluster label.
"""
max_cluster_label(m::DPM) = m.Q[1]

"""
    dp_mass(m::AbstractDPM)

Return the DP mass parameter.
"""
dp_mass(m::DPM) = m.α[1]

"""
    a0(m::AbstractDPM)

Return the shape parameter in the DP mass parameter's prior distribution.
"""
a0(m::DPM) = m.a0

"""
    b0(m::AbstractDPM)

Return the scale parameter in the DP mass parameter's prior distribution.
"""
b0(m::DPM) = m.b0

"""
    parent(m::AbstractDPM)

Return the parent DPM. 
"""
parent(m::DPM) = m

"""
    @forward_all_dpm_methods T.x

Extend all the methods originally defined on type `DPM` on type `T`, which call 
the relevant functions on the field `x`.

# One Example

```julia
struct CustomDPM
    gb::DPM
end
@forward_all_dpm_methods CustomDPM.gb
# Now dp_mass(m::Custom) call dp_mass(m.gb)
# Now cluster_labels(m::Custom) call cluster_labels(m.gb)
# etc.
```
"""
macro forward_dpm_methods(ex)
    return :(
        @forward $(ex) n_clusters, active_clusters, cluster_labels, 
        cluster_sizes, dp_mass, max_cluster_label, passive_clusters, parent
    )
end

# 3. Interface
# Any DPM specific block (sb) must implement these functions

"""
    logpredlik(sb::AbstractDPM_SB, gb::DPMGB, data, i, k)

Return ``\\log(y_i | y_{-i}, d_i = k, d_{-i})``. 
"""
function logpredlik(m::AbstractDPM, y, i, k)
    # Return log p(y[i] | y[-i], d[-i], d[i] = k)
    error("not implemented")
end

function update_hyperpars!(rng, m::AbstractDPM, y)
end

function update_suffstats!(m::AbstractDPM, y)
    # Update the sufficient statistics (if any) from scratch
    error("not implemented")
end

function update_suffstats!(m::AbstractDPM, y, i, k0, k1)
    # Update the sufficient statistics (if any) after `d[i]` changes from `k0` to `k1`
    error("not implemented")
end

# 4. Gibbs update logic

function update!(rng, m::AbstractDPM, y)
    update_hyperpars!(rng, m, y) # update the hyperparameters
    update_α!(rng, parent(m)) # update the DP mass parameter
    update_d!(rng, m, y) # update the cluster labels
end

function update_d!(rng, m::AbstractDPM, y)
    @unpack K, Q, A, P, d, n, τ, α = parent(m)
    update_suffstats!(m, y)
    for i in randperm!(rng, τ)
        d0 = d[i]
        d1 = first(P)
        p1 = log(α[1]) 
        p1 += logpredlik(m, y, i, d1)
        p1 -= log(-log(rand(rng)))
        for k in A
            p = logpredlik(m, y, i, k)
            p += log(n[k] - (d0 == k))
            p -= log(-log(rand(rng)))
            p > p1 && (d1 = k; p1 = p)
        end
        if d1 != d0
            (n[d0] -= 1) == 0 && (push!(P, d0); pop!(A, d0); K[1] -= 1)
            (n[d1] += 1) == 1 && (push!(A, d1); pop!(P, d1); K[1] += 1)
            isempty(P) && (push!(n, 0); push!(P, K[1] + 1); Q[1] += 1)
            update_suffstats!(m, y, i, d0, d1)
            d[i] = d1
        end
    end
    return nothing
end

function update_α!(rng, m::DPM)
    @unpack N, K, α, a0, b0 = m
    ϕ = rand(rng, Beta(α[1] + 1.0, N))
    ψ = 1.0 / (1.0 + N * (b0 - log(ϕ)) / (a0 + K[1] - 1.0))
    α[1] = rand(rng, Gamma(a0 + K[1] - (rand(rng) > ψ), 1.0 / (b0 - log(ϕ))))
    return nothing
end

include("normaldpm.jl")

end # module
