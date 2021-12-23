"""
    Skeleton(; N, M; 1_α0 = 2.0, b_α0 = 1.0)

Initialize a DPM Skeleton with `N` observations for training, `M` observations 
for prediction, and `Gamma(a_α0, 1 / b_α0)` prior distribution for the DP mass 
parameter.
"""
Base.@kwdef struct Skeleton
    # Data
    N::Int # sample size (training)
    M::Int # sample size (prediction)
    # Hyperparameters
    a_α0::Float64 = 1.0 # shape parameter in α's prior
    b_α0::Float64 = 1.0 # rate parameter in α's prior
    # Parameters
    α::Base.RefValue{Float64} = Ref(1.0) # dp mass parameter
    τ::Vector{Int} = collect(1:N)        # sample permutation
    d::Vector{Int} = ones(N)             # cluster labels
    # Transformed parameters
    f::Vector{Float64}  = zeros(M) # mixture density
    K::Base.RefValue{Int} = Ref(1) # number of clusters
    P::Set{Int} = Set(2)           # passive clusters
    A::Set{Int} = Set(1)           # active clusters
    n::Vector{Int} = [N]           # cluster sizes
end

"n_clusters(s::Skeleton) - Return the number of clusters."
n_clusters(s::Skeleton) = s.K[]

"cluster_sizes(s::Skeleton) - Return the size of each cluster."
component_sizes(s::Skeleton) = s.n

"active_clusters(s::Skeleton) - Return the set of active clusters."
active_clusters(s::Skeleton) = s.A

"passive_clusters(s::Skeleton) - Return the set of passive clusters."
passive_clusters(s::Skeleton) = s.P

"cluster_labels(s::Skeleton) - Return the cluster label of each sample point."
cluster_labels(s::Skeleton) = s.d

"predictive_pdf(s::Skeleton) - Return the predictive pdf at each grid point."
predictive_pdf(s::Skeleton) = s.f