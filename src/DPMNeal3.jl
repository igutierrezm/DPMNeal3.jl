module DPMNeal3

abstract type AbstractDPM end

using Distributions: Beta, Gamma 
using Parameters: @unpack
using Random: AbstractRNG, randperm, randperm!
export 
    # Types
    AbstractDPM, DPM,
    # Methods
    update!, cluster_labels, cluster_sizes, n_clusters, dp_mass, 
    active_clusters, passive_clusters, max_cluster_label

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
    function DPM(rng::AbstractRNG, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)
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
cluster_labels(m::AbstractDPM) = parent(m).d
cluster_labels(m::DPM) = m.d

"""
    cluster_sizes(m::AbstractDPM)

Return the cluster sizes.
"""
cluster_sizes(m::AbstractDPM) = parent(m).n
cluster_sizes(m::DPM) = m.n

"""
    active_clusters(m::AbstractDPM)

Return the set of active clusters.
"""
active_clusters(m::AbstractDPM) = parent(m).A
active_clusters(m::DPM) = m.A

"""
    passive_clusters(m::AbstractDPM)

Return a subset of the passive clusters.
"""
passive_clusters(m::AbstractDPM) = parent(m).P
passive_clusters(m::DPM) = m.P

"""
    n_clusters(m::AbstractDPM)

Return the number of active clusters.
"""
n_clusters(m::AbstractDPM) = parent(m).K[1]
n_clusters(m::DPM) = m.K[1]

"""
    max_cluster_label(m::AbstractDPM)

Return the largest (active) cluster label.
"""
max_cluster_label(m::AbstractDPM) = parent(m).Q[1]
max_cluster_label(m::DPM) = m.Q[1]

"""
    dp_mass(m::AbstractDPM)

Return the DP mass parameter.
"""
dp_mass(m::AbstractDPM) = parent(m).α[1]
dp_mass(m::DPM) = m.α[1]

# 3. Interface
# Any DPM specific block (sb) must implement these functions

"""
    parent(m::AbstractDPM)

Return the parent DPM.
"""
function parent(m::AbstractDPM)
    error("not implemented")
end

"""
    logpredlik(sb::AbstractDPM_SB, gb::DPMGB, data, i, k)

Return ``\\log(y_i | y_{-i}, x, d_i = k, d_{-i})``. 
"""
function logpredlik(m::AbstractDPM, data, i, k)
    # Return log p(y[i] | y[-i], d[-i], d[i] = k)
    error("not implemented")
end

function update_hyperpars!(rng::AbstractRNG, m::AbstractDPM, data)
end

function update_suffstats!(m::AbstractDPM, data)
    # Update the sufficient statistics (if any) from scratch
    error("not implemented")
end

function update_suffstats!(m::AbstractDPM, data, i, k0, k1)
    # Update the sufficient statistics (if any) after `d[i]` changes from `k0` to `k1`
    error("not implemented")
end

# 4. Gibbs update logic

function update!(rng::AbstractRNG, m::AbstractDPM, data)
    update_hyperpars!(rng, m, data) # update the hyperparameters
    update_d!(rng, m, data) # update the cluster labels
    update_α!(rng, m) # update the DP mass parameter
end

function update_d!(rng::AbstractRNG, m::AbstractDPM, data)
    @unpack K, Q, A, P, d, n, τ, α = parent(m)
    update_suffstats!(m, data)
    for i in randperm!(rng, τ)
        d0 = d[i]
        d1 = first(P)
        p1 = log(α[1]) 
        p1 += logpredlik(m, data, i, d1)
        p1 -= log(-log(rand(rng)))
        for k in A
            p = logpredlik(m, data, i, k)
            p += log(n[k] - (d0 == k))
            p -= log(-log(rand(rng)))
            p > p1 && (d1 = k; p1 = p)
        end
        if d1 != d0
            (n[d0] -= 1) == 0 && (push!(P, d0); pop!(A, d0); K[1] -= 1)
            (n[d1] += 1) == 1 && (push!(A, d1); pop!(P, d1); K[1] += 1)
            isempty(P) && (push!(n, 0); push!(P, K[1] + 1); Q[1] += 1)
            update_suffstats!(m, data, i, d0, d1)
            d[i] = d1
        end
    end
    return nothing
end

function update_α!(rng::AbstractRNG, m::AbstractDPM)
    @unpack N, K, α, a0, b0 = parent(m)
    ϕ = rand(rng, Beta(α[1] + 1.0, N))
    ψ = 1.0 / (1.0 + N * (b0 - log(ϕ)) / (a0 + K[1] - 1.0))
    α[1] = rand(rng, Gamma(a0 + K[1] - (rand(rng) > ψ), 1.0 / (b0 - log(ϕ))))
    return nothing
end

end # module
