"""
    Skeleton(m::AbstractModel)

Return the model skeleton.
"""
function skeleton(m::AbstractModel)
    error("not implemented")
end

"""
    in_sample_logpredlik(m::AbstractModel, i::Int, k::Int)

Return the log-predictive likelihood at the `i`th in-sample unit, given a 
cluster label equal to `j`, and the other responses and cluster labels in the 
sample.
"""
function in_sample_logpredlik(m::AbstractModel, i::Int, k::Int)
    error("not implemented")
end

"""
    out_of_sample_logpredlik(m::AbstractModel, i::Int, k::Int)

Return the log-predictive likelihood at the `i`th out-of-sample unit, given a 
cluster label equal to `j`, and the responses and cluster labels in the sample.
"""
function out_of_sample_logpredlik(m::AbstractModel, i::Int, k::Int)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractModel)

Update the suffstats from scratch.
"""
function update_suffstats!(m::AbstractModel)
    error("not implemented")
end

"""
    update_suffstats!(m::AbstractModel, i::Int, k0::Int, k1::Int)

Update the suffstats after the `i`th cluster label moves from `k0` to `k1`.
"""
function update_suffstats!(m::AbstractModel, i::Int, k::Int, l::Int)
    error("not implemented")
end
