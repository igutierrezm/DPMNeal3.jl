# Any new DPM, say MyDPM <: AbstractDPM, must extend DPM by composition
# and implement the following methods (with m::AbstractDPM replaced by MyDPM)

"""
    parent(m::AbstractDPM)

Return the parent DPM of `m`.
"""
function parent(m::AbstractDPM)
    error("not implemented")
end

"""
    logpredlik(m::AbstractDPM, data, i::Int, k::Int)

Return ``\\log(y_i | y_{-i}, x, d_i = k, d_{-i})``, where y_i, x_i and d_i are
the response, the covariates (if any) and the cluster label associated with 
the `i`th observation. Both x_i and y_i should be present in the dataset 
`data`.
"""
function logpredlik(m::AbstractDPM, data, i::Int, k::Int)
    error("not implemented")
end

"""
    update_hyperpars!(rng::AbstractRNG, m::AbstractDPM, data)

Update the kernel hyperparameters in `m`, given the dataset `data`.
"""
function update_hyperpars!(rng::AbstractRNG, m::AbstractDPM, data)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractDPM, data, i::Int, k0::Int, k1::Int)

Update the sufficient statistics in `m` from scratch, given the dataset `data`.
"""
function update_suffstats!(m::AbstractDPM, data)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractDPM, data, i::Int, k0::Int, k1::Int)

Update the sufficient statistics in `m` after the cluster label of the `i`th 
observation changes from `k0` to `k1`, given the dataset `data`.
"""
function update_suffstats!(m::AbstractDPM, data, i::Int, k0::Int, k1::Int)
    error("not implemented")
end
