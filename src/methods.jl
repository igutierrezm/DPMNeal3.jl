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

"""
    cluster_labels(m::AbstractDPM)

Return the cluster labels.
"""
cluster_labels(m::AbstractDPM) = parent(m).d

"""
    cluster_sizes(m::AbstractDPM)

Return the cluster sizes.
"""
cluster_sizes(m::AbstractDPM) = parent(m).n

"""
    active_clusters(m::AbstractDPM)

Return the set of active clusters.
"""
active_clusters(m::AbstractDPM) = parent(m).A

"""
    passive_clusters(m::AbstractDPM)

Return a subset of the passive clusters.
"""
passive_clusters(m::AbstractDPM) = parent(m).P

"""
    n_clusters(m::AbstractDPM)

Return the number of active clusters.
"""
n_clusters(m::AbstractDPM) = parent(m).K[1]

"""
    max_cluster_label(m::AbstractDPM)

Return the largest (active) cluster label.
"""
max_cluster_label(m::AbstractDPM) = parent(m).Q[1]

"""
    dp_mass(m::AbstractDPM)

Return the DP mass parameter, α.
"""
dp_mass(m::AbstractDPM) = parent(m).α[1]
