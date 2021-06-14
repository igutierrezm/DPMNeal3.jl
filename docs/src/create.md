# Creating new DPM models

To implement a new DPM, the first step is to define a subtype of `AbstractDPM` that extends the type `DPM` using composition, e.g.

```julia
struct MyDPM <: AbstractDPM
    parent::DPM
    # ... any desired fields
end
```

In this way, all the standard components of the DPM (e.g. the vector of cluster labels) will be stored in `parent`. Please, read carefully the documentation about `DPM` and its constructor/accessors in order to avoid redundancies.

Once `MyDPM` is defined, the last step is to implement the following methods:

```julia
parent(m::MyDPM)
update_suffstats!(m::MyDPM, data)
update_suffstats!(m::MyDPM, data, i::Int, k0::Int, k1::Int)
update_hyperpars!(rng::AbstractRNG, m::MyDPM, data)
logpredlik(m::MyDPM, data, i::Int, k::Int)
```

see the documentation of each method for more details. Once again, we note that `data` can have any type, as long as each method is correctly implemented.