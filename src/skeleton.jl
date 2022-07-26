"""
    Skeleton(; 
        N::Int,
        α::Base.RefVakue{Float64} = Ref(1.0),
        τ::Vector{Int} = collect(1:N),
        d::Vector{Int} = ones(Int, N),
        a0::Float64 = 2.0, 
        b0::Float64 = 1.0,
        check_args = true
    )

Initialize the skeleton of a DPM model applied to a sample of size `N`. Users 
can specify the shape `a0` and rate `b0` in α's prior distribution, as well as 
the initial state of the model parameters: the mass parameter `α`, the vector 
of cluster labels `d`, and the sample permutation `τ`. If `check_args` is 
true, the constructor will also check the arguments.
"""
struct Skeleton
    # Data
    N::Int # sample size
    # Hyperparameters
    a0::Float64 # shape parameter in α's prior
    b0::Float64 # rate parameter in α's prior
    # Parameters
    α::Base.RefValue{Float64} # dp mass parameter
    τ::Vector{Int} # sample permutation
    d::Vector{Int} # cluster labels
    # Transformed parameters
    K::Base.RefValue{Int} # number of clusters
    P::Set{Int} # passive clusters
    A::Set{Int} # active clusters
    n::Vector{Int} # cluster sizes
    function Skeleton(;
        N::Int,
        α::Base.RefValue{Float64} = Ref(1.0),
        τ::Vector{Int} = collect(1:N),
        d::Vector{Int} = ones(Int, N),
        a0::Float64 = 2.0,
        b0::Float64 = 1.0,
        check_args = true
    )
        check_args && check_arguments(N, α, d, τ, a0, b0)
        K = Ref(maximum(d))
        P = Set(1 + K[])
        A = Set(1:K[])
        n = zeros(Int, K[])
        for i in 1:N
            n[d[i]] += 1
        end
        new(N, a0, b0, α, τ, d, K, P, A, n)
    end
end

function check_arguments(N, α, d, τ, a0, b0)
    Nτ = length(τ)
    Nd = length(d)
    N != Nd && throw(DimensionMismatch("d must have lengths $N != $Nd"))
    N != Nτ && throw(DimensionMismatch("τ must have lengths $N != $Nτ"))
    N == 0 && throw(ArgumentError("d and τ cannot be empty vectors"))
    if sort(τ) != collect(1:Nτ)
        throw(DomainError("τ must be a permutation of 1:$Nd"))
    end
    if sort(unique(d)) != 1:maximum(d)
        throw(DomainError("d must have consecutive labels starting from 1"))
    end
    α[] <= 0 && throw(DomainError("α[] must be a positive value"))
    a0 <= 0 && throw(DomainError("a0 must be a positive value"))
    b0 <= 0 && throw(DomainError("b0 must be a positive value"))
end

"""
    n_clusters(s::Skeleton)

Return the number of clusters.
"""
function n_clusters(s::Skeleton)
    s.K[]
end

"""
    cluster_sizes(s::Skeleton; make_copy::Bool = true)

Return a make_copy of the size of each cluster. 
If `make_copy == false`, a binding to the current state of the 
variable is returned, avoiding unnecessary copies (use with caution).
"""
function cluster_sizes(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.n) : s.n
end

"""
    "active_clusters(s::Skeleton; make_copy::Bool = true)

Return a make_copy of the set of active clusters.
If `make_copy == false`, a binding to the current state of the 
variable is returned, avoiding unnecessary copies (use with caution).
"""
function active_clusters(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.A) : s.A
end

""" 
    passive_clusters(s::Skeleton; make_copy::Bool = true)
    
Return a make_copy of the set of passive clusters.
If `make_copy == false`, a binding to the current state of the 
variable is returned, avoiding unnecessary copies (use with caution).
"""
function passive_clusters(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.P) : s.P
end

"""
    cluster_labels(s::Skeleton; make_copy::Bool = true)
    
Return a make_copy of the cluster label of each observation. 
If `make_copy == false`, a binding to the current state of the 
variable is returned, avoiding unnecessary copies (use with caution).
"""
function cluster_labels(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.d) : s.d
end
