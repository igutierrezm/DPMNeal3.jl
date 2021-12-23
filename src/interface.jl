"""
    Skeleton(m::AbstractModel)

Return the model skeleton.
"""
function skeleton(m::AbstractModel)
    error("not implemented")
end

"""
    logpredlik(m::AbstractModel, i::Int, k::Int, prediction = false)

Return the log-predictive likelihood at the `i`th training outcome, given a 
cluster label equal to `j`, the other training units, and the other cluster 
labels. If `prediction == true`, the `i`th prediciton unit is used instead.
"""
function logpredlik(m::AbstractModel, i::Int, k::Int, prediction::Bool = false)
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
