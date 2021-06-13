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

Let us initialize a Normal DPM for a sample with `N = 1000` observations.

The first step is to set up the environment:

```julia
julia> using Random, DPMNeal3
julia> rng = MersenneTwister(1)
```

Then, we create the generic block (GB) of our DPM using `DPMGB()`:

```julia
julia> N = 1000
julia> gb = DPMGB(rng, N)
```

Finally, we create the specific block (GB) of our DPM using `DPM_SB_Normal()`:

```julia
julia> sb = DPM_SB_Normal()
```

## Performing a Gibbs update

In order to perform a Gibbs update, we need a training sample, e.g.

```julia
julia> y = randn(N)
julia> x = ones(N)
```

Given this sample, we can perform a Gibbs update as follows:

```julia
julia> update!(rng, sb, gb, y, x)
```

## Accesing the state of the chain

Once `m` have been updated, we can inspect the current state of the chain using the following accessors:

```julia
n_clusters(gb)      # return the number of active clusters
active_clusters(gb) # return the active clusters
cluster_labels(gb)  # return the cluster labels
cluster_sizes(gb)   # return the cluster sizes
dp_mass(gb)         # return the DP mass parameter
```
