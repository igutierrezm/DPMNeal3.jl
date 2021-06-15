```@meta
CurrentModule = DPMNeal3
```

# DPMNeal3.jl

[DPMNeal3.jl](https://github.com/igutierrezm/DPMNeal3.jl) provides an interface for Dirichlet process mixture (DPM) models in Julia using Neal's algorithm 3 [1]. Particularly, this package provides:

- A method for performing one iteration of Neal's algorithm 3.

- Several methods for accessing the current state of the chain (e.g. the vector of cluster labels).