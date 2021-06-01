using StatsBase

struct SpecificBlock
    N::Int # sample size
    G::Int # number of reduced groups
    γ::Vector{Vector{Float64}}
    z::Vector{Vector{Float64}}
    ng::Vector{Vector{Float64}}
    ybg::Vector{Vector{Float64}}
    s2g::Vector{Vector{Float64}}
    function SpecificBlock(N::Int, G::Int)
        γ = zeros(G)
        z = ones(N)
        m = [zeros(G)]
        yb = [zeros(G)]
        s2 = [zeros(G)]
        new(N, G, γ, z, m, yb, s2)
    end
end

# Compute the reduce group labels (g)
N = 10
D = 2
J = [3, 2]
x = [zeros(Int, D) for i in 1:N]
for i = 1:N, d = 1:D
    x[i][d] = rand(1:J[d])
end
g = StatsBase.competerank(x)

function update_suffstats!(sb::SpecificBlock, gb::GenericBlock, data)
    @unpack y, x = data
    @unpack N, K, n = gb
    @unpack m, yb, s2 = sb
    while length(yb):length(n)
        push!(m,  zeros(G))
        push!(yb, zeros(G))
        push!(s2, zeros(G))
    end
    for i = 1:N
        g = x[i]^γ[]
        n[x[]]
    end
end

function suffstats0!(s::SuffStats{A, B, C}, c::ChainState) where {A, B, C}
    @unpack m, y, x, K̄, N, J, n, r, ν, u, S = s
    @unpack r0, ν0, u0, S0 = m
    @unpack z, γ, O = c
    idx = Bool.(ones(N))
    @inbounds for j = 1:J, k ∈ 1:K̄
        n[k] = 0
        @. idx = (x^γ[x] == j) & (z == k)
        if any(idx)
            ysub = y[idx]
            njk = length(ysub)
            ym = mean(ysub)
            yv = cov(ysub, corrected = false)
            ν[j, k] = ν0 + njk
            r[j, k] = r0 + njk
            u[j, k] = (r0 * u0 + njk * ym) / r[j, k]
            S[j, k] = 
                (Matrix(S0) + reinterpret(Float64, njk * yv  + njk * r0 * (ym - u0) * (ym - u0)' / r[j, k])) |>
                x -> Symmetric(x) |>
                x -> cholesky(x)
        else
            ν[j, k] = ν0
            r[j, k] = r0
            u[j, k] .= u0
            S[j, k].factors .= S0.factors
        end
    end
    @inbounds @fastmath for i ∈ 1:N 
        n[z[i]] += 1 
    end
    return s
end