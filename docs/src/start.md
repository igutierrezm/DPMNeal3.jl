# Getting Started

## Installation

Install with the Julia package manager Pkg:

```julia
# Press ']' to enter the Pkg REPL mode.
pkg> add DPMNeal3
```
or
```julia
julia> using Pkg; 
julia> Pkg.add("DPMNeal3")
```

## Initializing a Normal DPM

Let us initialize a Normal DPM.

The first step is to set up the environment:

```julia
julia> using Random, DPMNeal3
julia> rng = MersenneTwister(1)
```

Then, we create a Normal DPM using `DPM()`:

```julia
julia> m = DPM(rng, NormalNIG())
```

## Performing a Gibbs update

In order to perform a Gibbs update, we need a training sample, e.g.

```julia
julia> data = randn(1000)
```

Given this sample, we can perform a Gibbs update as follows:

```julia
julia> update!(rng, m, data)
```

## Accesing the state of the chain

Once `m` have been updated, we can inspect the current state of the chain using the following accessors:

```julia
n_clusters(m)      # return the number of active clusters
active_clusters(m) # return the active clusters
cluster_labels(m)  # return the cluster labels
cluster_sizes(m)   # return the cluster sizes
dp_mass(m)         # return the DP mass parameter
```
