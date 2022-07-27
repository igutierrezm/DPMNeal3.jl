```@meta
CurrentModule = DPMNeal3
```

# DPMNeal3.jl

[DPMNeal3.jl](https://github.com/igutierrezm/DPMNeal3.jl) 
provides tools for the estimation of DPMs
[[1]](https://doi.org/10.1214/aos/1176346412) 
using Neal's algorithm 3 
[[2]](https://doi.org/10.1080/10618600.2000.10474879), including

1. A concrete type, `Skeleton`, to represent the *skeleton*[^1] of a DPM.
2. Functions for accessing the fields of a skeleton.
3. An abstract type, `AbstractDPM`, to represent a general DPM.
4. A function for updating the mass parameter of a general DPM as 
    in [[3]](https://doi.org/10.1080/01621459.1995.10476550).
5. A function for updating the vector of cluster indicators 
    of a general DPM using Neal's algorithm 3.

To use these tools with a new type `MyDPM <: AbstractDPM`, 
the type must specialize the methods in [Interface](@ref).

The `Skeleton` type and its accesors are described in
[The Skeleton type](@ref).

The `AbstractDPM` type and its mutators are described in
[The AbstractDPM type](@ref).

[^1]: Here, we use the term *skeleton* to denote all the generic elements of a
    DPM, such as the vector of cluster indicators.
