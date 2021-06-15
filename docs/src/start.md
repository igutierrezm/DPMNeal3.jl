# Getting Started

## Using a custom DPM

Let `m` be an object of type `MyDPM`, where `MyDPM <: AbstractDPM` is a datatype that conforms with the interface defined in this module. Then, `m` contains the current state of the chain associated with the Gibbs sampler described in Neal's algorithm 3. We can access the contents of `m` using the following accessors:

```julia
julia> DPMNeal3.n_clusters(m)      # return the number of active clusters
julia> DPMNeal3.active_clusters(m) # return the active clusters
julia> DPMNeal3.cluster_labels(m)  # return the cluster labels
julia> DPMNeal3.cluster_sizes(m)   # return the cluster sizes
julia> DPMNeal3.dp_mass(m)         # return the DP mass parameter
```

Now, let `data` be an object of type `MyData`, where `MyData <: Any` is the datatype of the sample expected by `MyDPM`. Then, we can perform one Gibbs update using `update!()`:

```julia
julia> rng = Random.MersenneTwister(1) # or any AbstractRNG object
julia> update!(rng, m, data)
```

Note that `data` can have any type, provided that `MyDPM` conforms with the interface.
