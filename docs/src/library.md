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
parent_dpm(::AbstractModel)
update_suffstats!(::AbstractModel, ::Any)
update_suffstats!(::AbstractModel, ::Any, ::Int, ::Int, ::Int)
update_hyperpars!
logpredlik
```