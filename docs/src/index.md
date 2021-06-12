```@meta
CurrentModule = DPMNeal3
```

# DPMNeal3

[DPMNeal3](https://github.com/igutierrezm/DPMNeal3.jl) provides an interface 
for Dirichlet process mixture (DPM) models in Julia using Neal's algorithm 3 
[1].

# Getting Started

## Installation

Install with the Julia package manager Pkg, just like any other registered Julia package:

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

Let us initialize a Normal DPM with `N = 1000` observations.

The first step is to set up the environment:

```julia
julia> using Random, DPMNeal3
julia> rng = MersenneTwister(1)
```

Then, we create the generic and specific blocks of this DPM `GenericBlock` and `SpecificBlock`, respectively:

```julia
julia> gb = GenericBlock(rng, 1000)
julia> sb = SpecificBlock(NormalDPM())
```

The generic block contains the components that are common to all DPMs (e.g. the cluster labels), while the specific block contains the components that are specific to Normal DPMs (e.g., the NIG kernel hyperparameters).

## Performing one Gibbs update

Now, suppose that our sample is contained in the random vector `data`. Then we can performe one Gibbs update as follows:
```julia
julia> update_chainstate!(rng, sb, gb, data)
```


Consider a DPM o in its so-called *stick-breaking* representation:

```math
\begin{aligned}
    y_i | \bm{d}, \bm{\theta}
    &\stackrel{\text{\tiny \it ind}}{\sim}
    q(\cdot | \theta_{d_i})
    \\
    d_i | \bm{w}
    &\stackrel{\text{\tiny \it iid}}{\sim}
    \text{Categorical}(\bm{w}),
    &
    \forall i
    &\in
    \mathcal{N} \equiv \{1, \ldots, N\}
    \\
    w_j
    &= 
    v_j \prod\nolimits_{z < j} (1 - v_z),
    &
    \forall j
    &\in
    \mathbb{N}
    \\
    v_j
    &\stackrel{\text{\tiny \it iid}}{\sim}
    \text{Beta}(1, \alpha)
    \\
    \theta_j
    &\stackrel{\text{\tiny \it iid}}{\sim}
    g_0(\cdot)
\end{aligned}
```


```@index
```

```@autodocs
Modules = [DPMNeal3]
```

# push!(LOAD_PATH,"../src/")
# rm -rf build; julia make.jl