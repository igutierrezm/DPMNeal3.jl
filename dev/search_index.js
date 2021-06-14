var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"[1] Neal, R. (2000). Markov chain sampling methods for Dirichlet process mixture models. Journal of Computational and Graphical Statistics, 9(2). https://doi.org/10.1080/10618600.2000.10474879","category":"page"},{"location":"theory/","page":"-","title":"-","text":"Consider a DPM o in its so-called stick-breaking representation:","category":"page"},{"location":"theory/","page":"-","title":"-","text":"beginaligned\n    y_i  bmd bmtheta\n    stackreltexttiny it indsim\n    q(cdot  theta_d_i)\n    \n    d_i  bmw\n    stackreltexttiny it iidsim\n    textCategorical(bmw)\n    \n    forall i\n    in\n    mathcalN equiv 1 ldots N\n    \n    w_j\n    = \n    v_j prodnolimits_z  j (1 - v_z)\n    \n    forall j\n    in\n    mathbbN\n    \n    v_j\n    stackreltexttiny it iidsim\n    textBeta(1 alpha)\n    \n    theta_j\n    stackreltexttiny it iidsim\n    g_0(cdot)\nendaligned","category":"page"},{"location":"theory/","page":"-","title":"-","text":"","category":"page"},{"location":"theory/","page":"-","title":"-","text":"Modules = [DPMNeal3]","category":"page"},{"location":"theory/#DPMNeal3.DPM","page":"-","title":"DPMNeal3.DPM","text":"DPM(rng::AbstractRNG, N; K0 = 1, a0 = 2.0, b0 = 4.0)\n\nInitialize a generic DPM with N observations, K0 initial clusters and a  Gamma(a0, b0) prior distribution for the DP mass parameter.\n\n\n\n\n\n","category":"type"},{"location":"theory/#DPMNeal3.active_clusters-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.active_clusters","text":"active_clusters(m::AbstractDPM)\n\nReturn the set of active clusters.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.cluster_labels-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.cluster_labels","text":"cluster_labels(m::AbstractDPM)\n\nReturn the cluster labels.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.cluster_sizes-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.cluster_sizes","text":"cluster_sizes(m::AbstractDPM)\n\nReturn the cluster sizes.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.dp_mass-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.dp_mass","text":"dp_mass(m::AbstractDPM)\n\nReturn the DP mass parameter, α.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.logpredlik-Tuple{AbstractDPM, Any, Int64, Int64}","page":"-","title":"DPMNeal3.logpredlik","text":"logpredlik(m::AbstractDPM, data, i::Int, k::Int)\n\nReturn log(y_i  y_-i x d_i = k d_-i), where yi, xi and di are the response, the covariates (if any) and the cluster label associated with  the ith observation. Both xi and y_i should be present in the dataset  data.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.max_cluster_label-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.max_cluster_label","text":"max_cluster_label(m::AbstractDPM)\n\nReturn the largest (active) cluster label.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.n_clusters-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.n_clusters","text":"n_clusters(m::AbstractDPM)\n\nReturn the number of active clusters.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.parent-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.parent","text":"parent(m::AbstractDPM)\n\nReturn the parent DPM of m.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.passive_clusters-Tuple{AbstractDPM}","page":"-","title":"DPMNeal3.passive_clusters","text":"passive_clusters(m::AbstractDPM)\n\nReturn a subset of the passive clusters.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.update_hyperpars!-Tuple{Random.AbstractRNG, AbstractDPM, Any}","page":"-","title":"DPMNeal3.update_hyperpars!","text":"update_hyperpars!(rng::AbstractRNG, m::AbstractDPM, data)\n\nUpdate the kernel hyperparameters in m, given the dataset data.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.update_suffstats!-Tuple{AbstractDPM, Any, Int64, Int64, Int64}","page":"-","title":"DPMNeal3.update_suffstats!","text":"update_suffstats!(m::AbstractDPM, data, i::Int, k0::Int, k1::Int)\n\nUpdate the sufficient statistics in m after the cluster label of the ith  observation changes from k0 to k1, given the dataset data.\n\n\n\n\n\n","category":"method"},{"location":"theory/#DPMNeal3.update_suffstats!-Tuple{AbstractDPM, Any}","page":"-","title":"DPMNeal3.update_suffstats!","text":"update_suffstats!(m::AbstractDPM, data, i::Int, k0::Int, k1::Int)\n\nUpdate the sufficient statistics in m from scratch, given the dataset data.\n\n\n\n\n\n","category":"method"},{"location":"theory/","page":"-","title":"-","text":"<!– push!(LOAD_PATH,\"../src/\") rm -rf build; julia make.jl –>","category":"page"},{"location":"create/#Creating-new-DPM-models","page":"Creating new DPMs","title":"Creating new DPM models","text":"","category":"section"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"To implement a new DPM, the first step is to define a subtype of AbstractDPM that extends the type DPM using composition, e.g.","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"struct MyDPM <: AbstractDPM\n    parent::DPM\n    # ... any desired fields\nend","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"In this way, all the standard components of the DPM (e.g. the vector of cluster labels) will be stored in parent. Please, read carefully the documentation about DPM and its constructor/accessors in order to avoid redundancies.","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"Once MyDPM is defined, the last step is to implement the following methods:","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"parent(m::MyDPM)\nupdate_suffstats!(m::MyDPM, data)\nupdate_suffstats!(m::MyDPM, data, i::Int, k0::Int, k1::Int)\nupdate_hyperpars!(rng::AbstractRNG, m::MyDPM, data)\nlogpredlik(m::MyDPM, data, i::Int, k::Int)","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"see the documentation of each method for more details. Once again, we note that data can have any type, as long as each method is correctly implemented.","category":"page"},{"location":"start/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"start/#Using-a-custom-DPM","page":"Getting Started","title":"Using a custom DPM","text":"","category":"section"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Let MyDPM <: AbstractDPM be a datatype that conforms with the interface defined in this module, and let m an object of type MyDPM. Then, m contains the current state of the chain associated with the Gibbs sampler described in Neal's algorithm 3. Then, we can access the contents of m using the following accessors:","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"julia> DPMNeal3.n_clusters(m)      # return the number of active clusters\njulia> DPMNeal3.active_clusters(m) # return the active clusters\njulia> DPMNeal3.cluster_labels(m)  # return the cluster labels\njulia> DPMNeal3.cluster_sizes(m)   # return the cluster sizes\njulia> DPMNeal3.dp_mass(m)         # return the DP mass parameter","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Now, let MyData <: Any be the datatype of the sample expected by MyDPM, and let data an object of type MyData. Then, we can perform one Gibbs update using update!():","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"julia> rng = Random.MersenneTwister(1) # or any AbstractRNG object\njulia> update!(rng, m, data)","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Note that data can have any type, provided that MyDPM conforms with the interface.","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Install with the Julia package manager Pkg:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"# Press ']' to enter the Pkg REPL mode.\npkg> add DPMNeal3","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"julia> using Pkg; \njulia> Pkg.add(\"DPMNeal3\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DPMNeal3","category":"page"},{"location":"#DPMNeal3","page":"Home","title":"DPMNeal3","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DPMNeal3 provides an interface for Dirichlet process mixture (DPM) models in Julia using Neal's algorithm 3 [1]. Particularly, this package provides:","category":"page"},{"location":"","page":"Home","title":"Home","text":"A method for performing one iteration of Neal's algorithm 3.\nSeveral methods for accessing the current state of the chain (e.g. the vector of cluster labels).","category":"page"}]
}
