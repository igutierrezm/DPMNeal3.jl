module DPMNeal3

using Distributions # Beta
using Random # randperm!

# All DPM models are subtypes of
abstract type AbstractModel end

# All DPM models must extend (by composition) the Skeleton type, defined in
include("skeleton.jl")

# All subtypes of AbstractModel must implement the interface described in
include("interface.jl")

# # The methods available for any AbstractModel are described in
# include("methods.jl")

# # The specific models are implemented in
# include("dpmnormal.jl")

end # module
