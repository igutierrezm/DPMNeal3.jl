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

Then, we create Normal DPM with `N = 1000` observations using `NormalDPM()`:

```julia
julia> gb = NormalDPM(rng; N = 1000)
```

## Performing a Gibbs update

In order to perform a Gibbs update, we need a training sample, e.g.

```julia
julia> y = randn(N)
julia> x = ones(N)
```

Given this sample, we can perform a Gibbs update using `gibbs_update!()`:

```julia
julia> gibbs_update!(rng, m, y, x)
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
