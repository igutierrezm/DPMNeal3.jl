using StatsBase
using SpecialFunctions

struct SpecificBlock
    G::Int
    v0::Float64
    r0::Float64
    u0::Float64
    s0::Float64
    v1::Vector{Vector{Float64}}
    r1::Vector{Vector{Float64}}
    u1::Vector{Vector{Float64}}
    s1::Vector{Vector{Float64}}
    γ::Vector{Float64}
    function SpecificBlock(G; v0 = 2, r0 = 1, u0 = 0.0, s0 = 1.0)
        γ = zeros(G)
        v1 = [zeros(G)]
        r1 = [zeros(G)]
        u1 = [zeros(G)]
        s1 = [zeros(G)]
        new(G, v0, r0, u0, s0, v1, r1, u1, s1, γ)
    end
end

function resize!(sb::SpecificBlock, n::Integer)
    @unpack G, v1, r1, u1, s1 = sb
    while length(v1) < n
        push!(v1, zeros(G))
        push!(r1, zeros(G))
        push!(u1, zeros(G))
        push!(s1, zeros(G))
    end    
end

function update_sb!(sb::SpecificBlock, gb::GenericBlock, data)
    @unpack y, x = data
    @unpack N, A, d, n = gb
    @unpack v1, r1, u1, s1, v0, r0, u0, s0, γ = sb
    length(v1) < length(n) && resize!(sb, length(n))
    for k in A
        v1[k] .= v0
        r1[k] .= r0
        u1[k] .= u0
        s1[k] .= s0
    end
    for i = 1:N
        zi = iszero(γ[x[i]]) ? x[i] : 1
        di = d[i]
        v1[di][zi] += 1
        rm = r1[di][zi] += 1
        um = u1[di][zi] = ((rm - 1) * u1[di][zi] + y[i]) / rm
        s1[di][zi] += (rm / (rm - 1)) * (y[i] - um)^2
    end
end

function update_sb!(sb::SpecificBlock, gb::GenericBlock, data, i, k1, k2)
    @unpack n = gb
    @unpack y, x = data
    @unpack v1, r1, u1, s1, γ = sb
    length(v1) < length(n) && resize!(sb, length(n))
    zi = iszero(γ[x[i]]) ? x[i] : 1

    # Modify cluster/group di/k2
    v1[k2][zi] += 1
    rm = r1[k2][zi] += 1
    um = u1[k2][zi] = ((rm - 1) * u1[k2][zi] + y[i]) / rm
    s1[k2][zi] += (rm / (rm - 1)) * (y[i] - um)^2

    # Modify cluster/group di/k1
    rm = r1[k1][zi]
    um = u1[k1][zi]
    s1[k1][zi] -= (rm / (rm - 1)) * (y[i] - um)^2
    u1[k1][zi]  = (rm * u1[k1][zi] - y[i]) / (rm - 1)
    v1[k1][zi] -= 1
    r1[k1][zi] -= 1
end

function log_pl(sb::SpecificBlock, gb::GenericBlock, data, i, k)
    @unpack y, x = data
    @unpack v1, r1, u1, s1, γ = sb
    @unpack d = gb
    yi = y[i]
    di = d[i]
    j = iszero(γ[x[i]]) ? x[i] : 1

    if di == k
        v̄1 = v1[k][j]
        r̄1 = r1[k][j]
        ū1 = u1[k][j]
        s̄1 = s1[k][j]
        v̄0 = v̄1 - 1
        r̄0 = r̄1 - 1
        ū0 = (r̄1 * ū1 - yi) / r̄0
        s̄0 = s̄1 - (r̄1 / r̄0) * (yi - ū1)^2
    else
        v̄0 = v1[k][j]
        r̄0 = r1[k][j]
        ū0 = u1[k][j]
        s̄0 = s1[k][j]
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
