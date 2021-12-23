function sample!(m::AbstractModel; iter = 4000, warmup = 2000, thin = 1)
    (; M, f) = skeleton(m)
    fchain = [zeros(M) for _ in 1:(iter - warmup) ÷ thin]
    for _ in 1:iter
        update!(m)
        fchain .= f
    end
    return fchain
end

function update!(m::AbstractModel)
    update_d!(m) # update the cluster labels
    update_α!(m) # update the mass parameter
    update_f!(m) # update the (conditional) posterior predictive density
end

function update_α!(m::AbstractModel)
    (; N, K, α, a_α0, b_α0) = skeleton(m)
    ϕ = rand(Beta(α[] + 1.0, N))
    ψ = 1.0 / (1.0 + N * (b_α0 - log(ϕ)) / (a_α0 + K[] - 1.0))
    α[] = rand(Gamma(a0 + K[] - (rand() > ψ), 1.0 / (b0 - log(ϕ))))
    return nothing
end

function update_d!(m::AbstractModel)
    (; α, τ, d, K, P, A, n) = skeleton(m)
    update_suffstats!(m)
    for i in randperm!(τ)
        k0 = d[i]
        k1 = first(P)
        p1 = logpredlik(m, i, k1) + log(α[]) - log(-log(rand()))
        for k in A
            p = logpredlik(m, i, k) + log(n[k] - (k == k0)) - log(-log(rand()))
            p > p1 && (k1 = k; p1 = p)
        end
        if k1 != k0
            (n[k0] -= 1) == 0 && (push!(P, k0); pop!(A, k0); K[] -= 1)
            (n[k1] += 1) == 1 && (push!(A, k1); pop!(P, k1); K[] += 1)
            isempty(P) && (push!(n, 0); push!(P, K[] + 1))
            update_suffstats!(m, i, k0, k1)
            d[i] = k1
        end
    end
    return nothing
end

function update_f!(m::AbstractModel)
    (; N, M, α, f, P, A, n) = skeleton(m)
    for i in 1:M
        f[i] = α[] * exp(logpredlik(m, i, first(P), true))
        for k in A 
            f[i] += n[k] * exp(logpredlik(m, i, k, true))
        end
        f[i] = 
        f[i] /= (N + α[])
    end
    return nothing
end
