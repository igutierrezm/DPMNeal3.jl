"""
    skeleton(m::AbstractDPM)::Skeleton

Return the DPM skeleton.
"""
function skeleton(m::AbstractDPM)::Skeleton
    error("not implemented")
end

"""
    logpredlik(m::AbstractDPM, i::Int, k::Int)

Return the log-predictive likelihood of the `i`th observation at its current
value, given that its cluster indicator is equal to `k`, and given the other 
responses and cluster indicators in the sample.
"""
function logpredlik(m::AbstractDPM, i::Int, k::Int)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractDPM)::Nothing

Update the sufficient statistics of each component from scratch.
"""
function update_suffstats!(m::AbstractDPM)::Nothing
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractDPM, i::Int, k0::Int, k1::Int)::Nothing

Update the sufficient statistics of the cluster `k0` and `k1` 
after the `i`th cluster indicator moves from `k0` to `k1`.
"""
function update_suffstats!(m::AbstractDPM, i::Int, k0::Int, k1::Int)::Nothing
    error("not implemented")
end
