"""
    Skeleton(; 
        N::Int,
        α::Base.RefVakue{Float64} = Ref(1.0),
        τ::Vector{Int} = collect(1:N),
        d::Vector{Int} = ones(Int, N),
        a0::Float64 = 2.0, 
        b0::Float64 = 1.0,
        check_args::Bool = true
    )

Initialize the skeleton of a DPM model of sample size `N`. Users can specify 
the shape `a0` and rate `b0` in α's prior distribution, as well as the initial 
state of the model parameters: the mass parameter `α`, the vector of cluster 
labels `d`, and the sample permutation `τ`. If `check_args`, the 
constructor will also check the arguments.
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
    P::Set{Int} # passive components
    A::Set{Int} # active components
    n::Vector{Int} # component sizes
    function Skeleton(;
        N::Int,
        α::Base.RefValue{Float64} = Ref(1.0),
        τ::Vector{Int} = collect(1:N),
        d::Vector{Int} = ones(Int, N),
        a0::Float64 = 2.0,
        b0::Float64 = 1.0,
        check_args::Bool = true
    )
        check_args && check_arguments(N, α, d, τ, a0, b0)
        A = Set(d)
        K = Ref(length(A))
        C = maximum(d) + 1
        P = setdiff(Set(1:C), A)
        n = zeros(Int, C)
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
    if minimum(d) < 1
        throw(DomainError("d must have positive elements"))
    end
    α[] <= 0 && throw(DomainError("α must refer to a positive value"))
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
    n_components(s::Skeleton)

Return the number of components.
"""
function n_components(s::Skeleton)
    length(s.A) + length(s.P)
end

"""
    component_sizes(s::Skeleton; make_copy::Bool = true)

Return a copy of the vector of component sizes. 
If `!make_copy`, no copy is made.
"""
function component_sizes(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.n) : s.n
end

"""
    active_components(s::Skeleton; make_copy::Bool = true)

Return a copy of the set of active components.
If `!make_copy`, no copy is made.
"""
function active_components(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.A) : s.A
end

""" 
    passive_components(s::Skeleton; make_copy::Bool = true)
    
Return a copy of the set of passive components.
If `!make_copy`, no copy is made.
"""
function passive_components(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.P) : s.P
end

"""
    cluster_indicators(s::Skeleton; make_copy::Bool = true)
    
Return a copy of the vector of cluster indicators. 
If `!make_copy`, no copy is made.
"""
function cluster_indicators(s::Skeleton; make_copy::Bool = true)
    make_copy ? copy(s.d) : s.d
end

"""
    mass_parameter(s::Skeleton)
    
Return the DP mass parameter.
"""
function mass_parameter(s::Skeleton)
    s.α[]
end
