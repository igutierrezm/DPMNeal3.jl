# Any new DPM, say MyDPM <: AbstractDPM, must extend DPM by composition
# and implement the following methods (with m::AbstractDPM replaced by MyDPM)

"""
    parent_dpm(m::AbstractDPM)

Return the parent DPM.
"""
function parent_dpm(m::AbstractDPM)
    error("not implemented")
end

"""
    logpredlik(m::AbstractDPM, data, i::Int, k::Int)

Return the log-pdf of the `i`th response at its current value, given the other 
responses, the other cluster labels, and the own cluster label fixed at `k`. 
The responses should be present in `data`.
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
    update_suffstats!(m::AbstractDPM, data)

Update the sufficient statistics in `m` from scratch, given the dataset `data`.
"""
function update_suffstats!(m::AbstractDPM, data)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractDPM, data, i::Int, k::Int, l::Int)

Update the sufficient statistics in `m` after the `i`th cluster label changes 
from `k` to `l`, given the dataset `data`.
"""
function update_suffstats!(m::AbstractDPM, data, i::Int, k::Int, l::Int)
    error("not implemented")
end
