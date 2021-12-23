module DPMNeal3

using Distributions
using SpecialFunctions

# All DPM models are subtypes of AbstractModel:
abstract type AbstractModel end

# All subtypes of AbstractModel must extend (by composition) the Skeleton type:
include("skeleton.jl")

# All subtypes of AbstractModel must implement the following interface:
include("interface.jl")

# The methods available for any AbstractModel are described here:
include("methods.jl")

# export 
#     # Types
#     AbstractModel, DPM,
#     # Methods
#     update!, cluster_labels, cluster_sizes, n_clusters, dp_mass, 
#     active_clusters, passive_clusters, cluster_capacity,
#     # Interface
#     parent_dpm, update_suffstats!, update_hyperpars!, logpredlik

end # module
