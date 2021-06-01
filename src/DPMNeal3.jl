module DPMNeal3

# 1. Imports/Exports

using Distributions: Beta, Gamma 
using Parameters: @unpack
using Random: randperm, randperm!
export GenericBlock, update!

# 2. Types

struct GenericBlock
    N::Int          # sample size
    K::Ref{Int}     # number of clusters
    α::Ref{Float64} # DP mass parameter
    τ::Vector{Int}  # ordering of the observations
    d::Vector{Int}  # cluster labels
    n::Vector{Int}  # cluster sizes
    P::Vector{Int}  # passive clusters
    A::Set{Int}     # active clusters
    a0::Float64     # shape parameter in p(α)
    b0::Float64     # rate parameter in p(α)
    function GenericBlock(rng, N::Int; K0::Int = 8, a0 = 2.0, b0 = 4.0)
        @assert N >= K0
        @assert K0 >= 0
        @assert a0 >= 0
        @assert b0 >= 0
        K = Ref(K0)
        α = Ref(1.0)
        τ = randperm(rng, N)
        d = rand(rng, 1:K0, N); d[1:K0] = 1:K0;
        n = [sum(d .== k) for k in 1:(K0 + 1)]
        P = [K0 + 1]
        A = Set(1:K0)
        new(N, K, α, τ, d, n, P, A, a0, b0)
    end
end

# 3. Interface 
# Any DPM specific block (sb) must implement these functions

function predloglik(sb, gb::GenericBlock, data, i, k)
    # Return log p(y[i] | y[-i], d[-i], d[i] = k)
    error("not implemented")
end

function update_suffstats!(sb, gb::GenericBlock, data)
    # Recompute the sufficient statistics from scratch
    error("not implemented")
end

function update_suffstats!(sb, gb::GenericBlock, data, i, k0, k1)
    # Recompute the sufficient statistics after `d[i]` changes from `k0` to `k1`
    error("not implemented")
end

# 4. Gibbs update logic

function update!(rng, sb, gb::GenericBlock, data)
    update_d!(rng, sb, gb, data)
    update_α!(rng, gb)
end

function update_d!(rng, sb, gb::GenericBlock, data)
    @unpack A, P, d, n, τ, α = gb
    update_suffstats!(sb, gb, data)

    for i in randperm!(rng, τ)
        d0 = d[i]
        d1 = P[end]
        p1 = predloglik(sb, gb, data, i, d1) + log(α[]) - log(-log(rand(rng)))
        for k in A
            p = predloglik(sb, gb, data, i, k) 
            p += log(n[k] - (d0 == k))
            p -= log(-log(rand(rng)))
            p > p1 && (d1 = k; p1 = p)
        end
        if d1 != d0
            update_suffstats!(sb, gb, data, i, d0, d1)
            (n[d0] -= 1) == 0 && (push!(P, d0); pop!(A, d0); K[] -= 1)
            (n[d1] += 1) == 1 && (push!(A, d1); pop!(P, d1); K[] += 1)
            isempty(P) && (push!(n, 0); push!(P, K[] + 1))
            d[i] = d1
        end
    end
end

function update_α!(rng, gb::GenericBlock)
    @unpack N, K, α, a0, b0 = gb
    ϕ = rand(rng, Beta(α[] + 1.0, N))
    ψ = 1.0 / (1.0 + N * (b0 - log(ϕ)) / (a0 + M[] - 1.0))
    α[] = rand(rng, Gamma(a0 + K[] - (rand(rng) > ψ), 1.0 / (b0 - log(ϕ))))
end

end # module
