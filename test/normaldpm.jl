struct NormalDPM <: AbstractDPM
    parent::DPM
    v0::Float64
    r0::Float64
    u0::Float64
    s0::Float64
    v1::Vector{Float64}
    r1::Vector{Float64}
    u1::Vector{Float64}
    s1::Vector{Float64}
    function NormalDPM(
            rng::AbstractRNG, 
            N::Int; 
            K0::Int = 1, 
            a0::Float64 = 2.0, 
            b0::Float64 = 4.0, 
            v0::Float64 = 2.0, 
            r0::Float64 = 1.0, 
            u0::Float64 = 0.0, 
            s0::Float64 = 1.0
        )
        parent = DPM(rng, N; K0, a0, b0)
        v1 = [v0]
        r1 = [r0]
        u1 = [u0]
        s1 = [s0]
        new(parent, v0, r0, u0, s0, v1, r1, u1, s1)
    end
end

function parent_dpm(m::NormalDPM)
    m.parent
end

function update_hyperpars!(rng::AbstractRNG, m::NormalDPM, data)
end

function add_cluster!(m::NormalDPM)
    @unpack v0, r0, u0, s0, v1, r1, u1, s1 = m
    push!(v1, v0)
    push!(r1, r0)
    push!(u1, u0)
    push!(s1, s0)
end

function update_suffstats!(m::NormalDPM, data)
    N = length(data)
    d = cluster_labels(m)
    A = active_clusters(m)
    Q = cluster_capacity(m)
    @unpack v0, r0, u0, s0, v1, r1, u1, s1 = m
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
        u1[k] = ((r1[k] - 1) * u1[k] + data[i]) / r1[k]
        s1[k] += (r1[k] / (r1[k] - 1)) * (data[i] - u1[k])^2
    end
end

function update_suffstats!(m::NormalDPM, data, i::Int, k1::Int, k2::Int)
    Q = cluster_capacity(m)
    @unpack v1, r1, u1, s1 = m
    while length(v1) < Q
        add_cluster!(m)
    end

    # Modify suffstats for cluster k1
    s1[k1] -= (r1[k1] / (r1[k1] - 1)) * (data[i] - u1[k1])^2
    u1[k1]  = (r1[k1] * u1[k1] - data[i]) / (r1[k1] - 1)
    v1[k1] -= 1
    r1[k1] -= 1

    # Modify suffstats for cluster k2
    v1[k2] += 1
    r1[k2] += 1
    u1[k2] = ((r1[k2] - 1) * u1[k2] + data[i]) / r1[k2]
    s1[k2] += (r1[k2] / (r1[k2] - 1)) * (data[i] - u1[k2])^2
end

function logpredlik(m::NormalDPM, data, i::Int, k::Int)
    @unpack v1, r1, u1, s1 = m
    d = cluster_labels(m)
    yi = data[i]
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

    if iszero(r̄0)
        return - 0.5v̄1 * log(s̄1) + loggamma(v̄1 / 2) - 0.5 * log(π)
    else
        return (
            0.5v̄0 * log(s̄0) -
            0.5v̄1 * log(s̄1) +
            loggamma(v̄1 / 2) -
            loggamma(v̄0 / 2) +
            0.5 * log(r̄0 / r̄1) -
            0.5 * log(π)
        )        
    end
end
