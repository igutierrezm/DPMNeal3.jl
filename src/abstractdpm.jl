"""
    AbstractDPM

An abstract type representing a general DPM
"""
abstract type AbstractDPM end

"""
    update_mass_parameter!(m::AbstractDPM)::Nothing

Update the DP mass parameter following Escobar & West (1995).
"""
function update_mass_parameter!(m::AbstractDPM)::Nothing
    (; N, K, α, a0, b0) = skeleton(m)
    if N > 0
        ϕ = rand(Beta(α[] + 1.0, N))
        ψ = 1.0 / (1.0 + N * (b0 - log(ϕ)) / (a0 + K[] - 1.0))
        α[] = rand(Gamma(a0 + K[] - (rand() > ψ), 1.0 / (b0 - log(ϕ))))
    else
        α[] = rand(Gamma(a0, 1 / b0))
    end
    return nothing
end

"""
    update_cluster_indicators!(m::AbstractDPM)::Nothing

Update the DP mass parameter following Neal (2000) algorithm 3.
"""
function update_cluster_indicators!(m::AbstractDPM)::Nothing
    (; α, τ, d, K, P, A, n) = skeleton(m)
    update_suffstats!(m)
    for i in randperm!(τ)
        k0 = d[i]
        k1 = first(P)
        p1 = logpredlik(m, i, k1) +
            log(α[]) - log(-log(rand()))
        for k in A
            p = logpredlik(m, i, k) + 
                log(n[k] - (k == k0)) - log(-log(rand()))
            p > p1 && (k1 = k; p1 = p)
        end
        if k1 != k0
            k1 > length(n) && push!(n, 0)
            (n[k0] -= 1) == 0 && (push!(P, k0); pop!(A, k0); K[] -= 1)
            (n[k1] += 1) == 1 && (push!(A, k1); pop!(P, k1); K[] += 1)
            isempty(P) && (push!(n, 0); push!(P, K[] + 1))
            update_suffstats!(m, i, k0, k1)
            d[i] = k1
        end
    end
    return nothing
end
