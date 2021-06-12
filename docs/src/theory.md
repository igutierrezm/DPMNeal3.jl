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

<!-- push!(LOAD_PATH,"../src/")
rm -rf build; julia make.jl -->