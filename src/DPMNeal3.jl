module DPMNeal3

using Distributions: Beta, Gamma 
using Parameters: @unpack
using Random: randperm, randperm!
export GenericBlock, update!

"""
    GenericBlock(rng::AbstractRNG, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)

Initialize the generic block of a DPM, where `N` is the sample size, `K0` is 
the initial number of clusters, and (`a0`, `b0`) determine the prior 
distribution of the DP mass parameter, a `Gamma(a0, 1 / b0)` 
distribution. 
"""
struct GenericBlock
    N::Int          # sample size
    K::Ref{Int}     # number of clusters
    α::Ref{Float64} # DP mass parameter
    τ::Vector{Int}  # ordering of the observations
    d::Vector{Int}  # cluster labels
    n::Vector{Int}  # cluster sizes
    P::Set{Int}     # passive clusters
    A::Set{Int}     # active clusters
    a0::Float64     # shape parameter in α's prior
    b0::Float64     # rate parameter in α's prior
    function GenericBlock(rng, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)
        @assert N >= K0
        @assert K0 >= 0
        @assert a0 >= 0
        @assert b0 >= 0
        K = Ref(K0)
        α = Ref(1.0)
        τ = randperm(rng, N)
        d = rand(rng, 1:K0, N); d[1:K0] = 1:K0;
        n = [sum(d .== k) for k in 1:(K0 + 1)]
        P = Set(K0 + 1)
        A = Set(1:K0)
        new(N, K, α, τ, d, n, P, A, a0, b0)
    end
end

# 3. Interface 
# Any DPM specific block (sb) must implement these functions

"""
    logh(sb, gb::GenericBlock, data, i)

Return the log of ``h(y_i) := \\int q(y_i | \\theta_1) g(\\theta_1) d\\theta_1``
"""
function logh(sb, gb::GenericBlock, data, i)
    # Return the log of h(y[i]) := ∫ q(y[i] | θ[1]) g(θ[1]) dθ[1]
    error("not implemented")
end

function logq(sb, gb::GenericBlock, data, i, k)
    # Return log p(y[i] | y[-i], d[-i], d[i] = k)
    error("not implemented")
end

function update_sb!(sb, gb::GenericBlock, data)
    # Update sb from scratch
    error("not implemented")
end

function update_sb!(sb, gb::GenericBlock, data, i, k0, k1)
    # Update sb after `d[i]` changes from `k0` to `k1`
    error("not implemented")
end

# 4. Gibbs update logic

function update!(rng, sb, gb::GenericBlock, data)
    update_d!(rng, sb, gb, data) # update the cluster labels
    update_α!(rng, gb) # update α
end

function update_d!(rng, sb, gb::GenericBlock, data)
    @unpack K, A, P, d, n, τ, α = gb
    update_sb!(sb, gb, data)
    for i in randperm!(rng, τ)
        d0 = d[i]
        d1 = first(P)
        p1 = log(α[]) 
        p1 += logh(sb, gb, data, i)
        p1 -= log(-log(rand(rng)))
        for k in A
            p = logq(sb, gb, data, i, k)
            p += log(n[k] - (d0 == k))
            p -= log(-log(rand(rng)))
            p > p1 && (d1 = k; p1 = p)
        end
        if d1 != d0
            update_sb!(sb, gb, data, i, d0, d1)
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
    ψ = 1.0 / (1.0 + N * (b0 - log(ϕ)) / (a0 + K[] - 1.0))
    α[] = rand(rng, Gamma(a0 + K[] - (rand(rng) > ψ), 1.0 / (b0 - log(ϕ))))
end

end # module
