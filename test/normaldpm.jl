Base.@kwdef struct DPMNormal <: AbstractModel
    # Data
    M::Int = 50
    y::Vector{Float64}
    ygrid::Vector{Float64} = range(minimum(y), maximum(y), length = M)
    # Hyperparameters
    v0::Float64 = 1.0
    r0::Float64 = 1.0
    u0::Float64 = 1.0
    s0::Float64 = 1.0
    a_α0::Float64 = 1.0
    b_α0::Float64 = 1.0
    # Transformed parameters
    v1::Vector{Float64} = zeros(1)
    r1::Vector{Float64} = zeros(1)
    u1::Vector{Float64} = zeros(1)
    s1::Vector{Float64} = zeros(1)
    # Skeleton
    skl::Skeleton = Skeleton(N = length(y), M = length(ygrid); a_α0, b_α0)
end

function skeleton(m::DPMNormal)
    return m.skl
end

function out_of_sample_logpredlik(m::DPMNormal, i::Int, k::Int)
    (; N, ygrid, v1, r1, u1, s1) = m
    yi = ygrid[i]
    v̄0 = v1[k]
    r̄0 = r1[k]
    ū0 = u1[k]
    s̄0 = s1[k]
    v̄1 = v̄0 + 1
    r̄1 = r̄0 + 1
    ū1 = (r̄0 * ū0 + yi) / r̄1
    s̄1 = s̄0 + (r̄1 / r̄0) * (yi - ū1)^2
    return common_logpredlik(v̄0, r̄0, ū0, s̄0, v̄1, r̄1, ū1, s̄1)
end

function in_sample_logpredlik(m::DPMNormal, i::Int, k::Int)
    (; N, y, v1, r1, u1, s1, skl) = m
    d = cluster_labels(skl)
    yi = y[i]
    di = d[i]
    if di == k
        v̄1 = v1[k]
        r̄1 = r1[k]
        ū1 = u1[k]
        s̄1 = s1[k]
        v̄0 = v̄1 - 1
        r̄0 = r̄1 - 1
        ū0 = (r̄1 * ū1 - yi) / r̄0
        s̄0 = s̄1 - (r̄1 / r̄0) * (yi - ū1)^2
    else
        v̄0 = v1[k]
        r̄0 = r1[k]
        ū0 = u1[k]
        s̄0 = s1[k]
        v̄1 = v̄0 + 1
        r̄1 = r̄0 + 1
        ū1 = (r̄0 * ū0 + yi) / r̄1
        s̄1 = s̄0 + (r̄1 / r̄0) * (yi - ū1)^2
    end
    return common_logpredlik(v̄0, r̄0, ū0, s̄0, v̄1, r̄1, ū1, s̄1)
end

function common_logpredlik(v̄0, r̄0, ū0, s̄0, v̄1, r̄1, ū1, s̄1)
    return (
        0.5v̄0 * log(s̄0) -
        0.5v̄1 * log(s̄1) +
        loggamma(v̄1 / 2) -
        loggamma(v̄0 / 2) +
        0.5 * log(r̄0 / r̄1) -
        0.5 * log(π)
    )        
end

function update_suffstats!(m::DPMNormal)
    (; N, v0, r0, u0, s0, v1, r1, u1, s1, skl) = m
    d = cluster_labels(skl)
    A = active_clusters(skl)
    Q = first(passive_clusters(skl))
    while length(v1) < Q
        add_cluster!(m)
    end
    for k in A
        v1[k] = v0
        r1[k] = r0
        u1[k] = u0
        s1[k] = s0
    end
    for i = 1:N
        k = d[i]
        v1[k] += 1
        r1[k] += 1
        u1[k] = ((r1[k] - 1) * u1[k] + y[i]) / r1[k]
        s1[k] += (r1[k] / (r1[k] - 1)) * (y[i] - u1[k])^2
    end
end

function update_suffstats!(m::DPMNormal, i::Int, k0::Int, k1::Int)
    (; y, v1, r1, u1, s1, skl) = m
    Q = first(passive_clusters(skl))
    while length(v1) < K0
        add_cluster!(m)
    end

    # Modify suffstats for cluster k0
    s1[k0] -= (r1[k0] / (r1[k0] - 1)) * (y[i] - u1[k0])^2
    u1[k0]  = (r1[k0] * u1[k0] - y[i]) / (r1[k0] - 1)
    v1[k0] -= 1
    r1[k0] -= 1

    # Modify suffstats for cluster k1
    v1[k1] += 1
    r1[k1] += 1
    u1[k1] = ((r1[k1] - 1) * u1[k1] + y[i]) / r1[k1]
    s1[k1] += (r1[k1] / (r1[k1] - 1)) * (y[i] - u1[k1])^2
end

function add_cluster!(m::DPMNormal)
    (; v0, r0, u0, s0, v1, r1, u1, s1) = m
    push!(v1, v0)
    push!(r1, r0)
    push!(u1, u0)
    push!(s1, s0)
end
