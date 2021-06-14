"""
    DPM(rng::AbstractRNG, N; K0 = 1, a0 = 2.0, b0 = 4.0)

Initialize a generic DPM with `N` observations, `K0` initial clusters and a 
`Gamma(a0, b0)` prior distribution for the DP mass parameter.
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

# Accessors

cluster_labels(m::DPM) = m.d
cluster_sizes(m::DPM) = m.n
active_clusters(m::DPM) = m.A
passive_clusters(m::DPM) = m.P
n_clusters(m::DPM) = m.K[1]
max_cluster_label(m::DPM) = m.Q[1]
dp_mass(m::DPM) = m.α[1]
