module DPMNeal3

using Distributions # Beta
using Random # randperm!

# We use the abstract type `AbstractDPM` for represeting a general DPM.
# The `AbstractDPM` type and its mutators are defined here:
include("abstractdpm.jl")

# We use the concrete type `Skeleton` for represeting a DPM skeleton[^1].
# The `Skeleton` type, its constructors and accessors and implemented here:
include("skeleton.jl")

end # module

# [^1]: 
# Here, we use the term skeleton to denote all the generic elements of a DPM, 
# such as the vector of cluster indicators.