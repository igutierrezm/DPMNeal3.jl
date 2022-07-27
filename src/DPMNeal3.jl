module DPMNeal3

using Distributions # Beta
using Random # randperm!

# We use the (abstract) type `AbstractDPM` to represent a general DPM.
# The type and its mutators are defined here:
include("abstractdpm.jl")

# We use the (concrete) type `Skeleton` to represet a DPM skeleton[^1].
# The type, its constructors and accessors and implemented here:
include("skeleton.jl")

include("interface.jl")

end # module

# [^1]: 
# A combination of the standard elements of a DPM, such as the mass parameter.