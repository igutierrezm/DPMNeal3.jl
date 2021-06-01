module DPMNeal3

using Parameters: @with_kw, @unpack
using Random: AbstractRNG, randperm!

@with_kw struct GenericBlock
    N::Int                        # sample size
    α::Ref{Float64} = Ref(1.0)    # DP mass parameter
    τ::Vector{Int} = collect(1:N) # ordering of the observations
    d::Vector{Int} = ones(N)      # cluster labels
    n::Vector{Int} = [N]          # cluster sizes
    P::Vector{Int} = [2]          # passive clusters
    A::Set{Int} = [1]             # active clusters
    a0::Float64 = 2               # shape parameter in α's prior
    b0::Float64 = 4               # rate parameter in α's prior
end

abstract type AbstractSpecificBlock end
abstract type AbstractDPMModel end

struct DPMModel{SpecificBlock <: AbstractSpecificBlock} <: AbstractDPMModel
    sb::SpecificBlock
    gb::GenericBlock
end

logpdf(m::AbstractDPMModel, y, i, j) = error("Not implemented")
update_suffstats!(m::AbstractDPMModel, y) = error("not implemented")
update_suffstats!(m::AbstractDPMModel, y, i, j1, j2) = error("not implemented")

function update!(m::AbstractDPMModel, rng, y)
    # Compute the sufficient statistics from scratch
    @unpack A, P, d, n, τ, α = m.gb
    update_suffstats!(m, y)

    # Update each cluster label
    for i in randperm!(rng, τ)
        # Remove one observation from cluster d0
        d0 = d[i]
        n[d0] -= 1
        if n[d0] == 0
            pop!(A, d0)
            push!(P, d0)
        end

        # Update d[i]
        d1 = P[end]
        p1 = logpdf(m, y, i, d1) + log(α[]) - log(-log(rand(rng)))
        for k in A
            p = logpdf1(m, y, i, k) + log(n[k]) - log(-log(rand(rng)))
            p < p1 && continue
            d1 = k
            p1 = p
        end
        d[i] = d1

        # Resize n if necessary
        if length(n) < d1
            push!(n, 0)
        end

        # Add one observation to cluster d1
        if n[d1] == 0
            push!(A, d1)
            last_popped = pop!(P, d1)
            isempty(P) && push!(P, last_popped + 1)
        end
        n[d1] += 1

        # Update the sufficient statistics as d[i] changes
        update_suffstats!(m, y, i, d0, d1)
    end
end

end # module

# N(m::AbstractDPMModel) = m.gb.N
# d(m::AbstractDPMModel) = m.gb.d
# n(m::AbstractDPMModel) = m.gb.n
# D(m::AbstractDPMModel) = m.gb.D
# α(m::AbstractDPMModel) = m.gb.α

