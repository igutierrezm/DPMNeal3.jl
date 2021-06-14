module DPMNeal3

abstract type AbstractDPM end

using Distributions: Beta, Gamma 
using Parameters: @unpack
using Random: AbstractRNG, randperm, randperm!

include("dpm.jl")
include("interface.jl")
include("methods.jl")

export 
    # Types
    AbstractDPM, DPM,
    # Methods
    update!, cluster_labels, cluster_sizes, n_clusters, dp_mass, 
    active_clusters, passive_clusters, max_cluster_label

end # module
