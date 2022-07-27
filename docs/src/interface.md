# Interface

In order to use the tools offered by this package with a new DPM model of
type `MyDPM <: AbstractDPM`, the following methods must be specialized
to `MyDPM`:[^1]

```@docs
DPMNeal3.skeleton
DPMNeal3.logpredlik
DPMNeal3.update_suffstats!
```

[^1]: In a nutshell, we should explain how to compute a predictive 
    log-likelihood (given the cluster indicators) and how to 
    access the skeleton.