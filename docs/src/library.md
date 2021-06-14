# Library

## Types

```@docs
DPM
```

## Accessors

```@docs
n_clusters
active_clusters
passive_clusters
cluster_capacity
cluster_labels
cluster_sizes
dp_mass
```

## Interface

```@docs
parent_dpm(::AbstractDPM)
update_suffstats!(::AbstractDPM, ::Any)
update_suffstats!(::AbstractDPM, ::Any, ::Int, ::Int, ::Int)
update_hyperpars!
logpredlik
```