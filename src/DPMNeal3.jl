module DPMNeal3

using Distributions
using Random
using SpecialFunctions

# All DPM models are subtypes of AbstractModel:
abstract type AbstractModel end

# All subtypes of AbstractModel must extend (by composition) the Skeleton type:
include("skeleton.jl")

# All subtypes of AbstractModel must implement the following interface:
include("interface.jl")

# The methods available for any AbstractModel are described here:
include("methods.jl")

# The specific models are implemented here:
include("dpmnormal.jl")

end # module
