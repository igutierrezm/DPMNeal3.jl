var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"[1] Neal, R. (2000). Markov chain sampling methods for Dirichlet process mixture models. Journal of Computational and Graphical Statistics, 9(2). https://doi.org/10.1080/10618600.2000.10474879","category":"page"},{"location":"library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/#Types","page":"Library","title":"Types","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"DPM","category":"page"},{"location":"library/#DPMNeal3.DPM","page":"Library","title":"DPMNeal3.DPM","text":"DPM(rng::AbstractRNG, N; K0 = 1, a0 = 2.0, b0 = 4.0)\n\nInitialize a generic DPM with N observations, K0 initial clusters and a  Gamma(a0, b0) prior distribution for the DP mass parameter.\n\n\n\n\n\n","category":"type"},{"location":"library/#Accessors","page":"Library","title":"Accessors","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"n_clusters\nactive_clusters\npassive_clusters\ncluster_capacity\ncluster_labels\ncluster_sizes\ndp_mass","category":"page"},{"location":"library/#DPMNeal3.n_clusters","page":"Library","title":"DPMNeal3.n_clusters","text":"n_clusters(m::AbstractDPM)\n\nReturn the current number of active clusters.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.active_clusters","page":"Library","title":"DPMNeal3.active_clusters","text":"active_clusters(m::AbstractDPM)\n\nReturn the current set of active clusters.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.passive_clusters","page":"Library","title":"DPMNeal3.passive_clusters","text":"passive_clusters(m::AbstractDPM)\n\nReturn (a subset of) the current set of passive clusters.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.cluster_capacity","page":"Library","title":"DPMNeal3.cluster_capacity","text":"cluster_capacity(m::AbstractDPM)\n\nReturn the current cluster storage capacity.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.cluster_labels","page":"Library","title":"DPMNeal3.cluster_labels","text":"cluster_labels(m::AbstractDPM)\n\nReturn the current cluster labels.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.cluster_sizes","page":"Library","title":"DPMNeal3.cluster_sizes","text":"cluster_sizes(m::AbstractDPM)\n\nReturn the current cluster sizes.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.dp_mass","page":"Library","title":"DPMNeal3.dp_mass","text":"dp_mass(m::AbstractDPM)\n\nReturn the current DP mass parameter.\n\n\n\n\n\n","category":"function"},{"location":"library/#Interface","page":"Library","title":"Interface","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"parent_dpm(::AbstractDPM)\nupdate_suffstats!(::AbstractDPM, ::Any)\nupdate_suffstats!(::AbstractDPM, ::Any, ::Int, ::Int, ::Int)\nupdate_hyperpars!\nlogpredlik","category":"page"},{"location":"library/#DPMNeal3.parent_dpm-Tuple{AbstractDPM}","page":"Library","title":"DPMNeal3.parent_dpm","text":"parent_dpm(m::AbstractDPM)\n\nReturn the parent DPM.\n\n\n\n\n\n","category":"method"},{"location":"library/#DPMNeal3.update_suffstats!-Tuple{AbstractDPM, Any}","page":"Library","title":"DPMNeal3.update_suffstats!","text":"update_suffstats!(m::AbstractDPM, data)\n\nUpdate the sufficient statistics in m from scratch, given the dataset data.\n\n\n\n\n\n","category":"method"},{"location":"library/#DPMNeal3.update_suffstats!-Tuple{AbstractDPM, Any, Int64, Int64, Int64}","page":"Library","title":"DPMNeal3.update_suffstats!","text":"update_suffstats!(m::AbstractDPM, data, i::Int, k::Int, l::Int)\n\nUpdate the sufficient statistics in m after the ith cluster label changes  from k to l, given the dataset data.\n\n\n\n\n\n","category":"method"},{"location":"library/#DPMNeal3.update_hyperpars!","page":"Library","title":"DPMNeal3.update_hyperpars!","text":"update_hyperpars!(rng::AbstractRNG, m::AbstractDPM, data)\n\nUpdate the kernel hyperparameters in m, given the dataset data.\n\n\n\n\n\n","category":"function"},{"location":"library/#DPMNeal3.logpredlik","page":"Library","title":"DPMNeal3.logpredlik","text":"logpredlik(m::AbstractDPM, data, i::Int, k::Int)\n\nReturn the log-pdf of the ith response at its current value, given the other  responses, the other cluster labels, and the own cluster label fixed at k.  The responses should be present in data.\n\n\n\n\n\n","category":"function"},{"location":"create/#Creating-new-DPM-models","page":"Creating new DPMs","title":"Creating new DPM models","text":"","category":"section"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"To implement a new DPM, the first step is to define a subtype of AbstractDPM that extends the type DPM using composition, e.g.","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"struct MyDPM <: AbstractDPM\n    parent::DPM\n    # ... any desired fields\nend","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"In this way, all the standard components of the DPM (e.g. the vector of cluster labels) will be stored in parent. Please, read carefully the documentation about DPM and its constructor/accessors in order to avoid redundancies.","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"Once MyDPM is defined, the last step is to specialize the following methods to m::MyDPM:","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"parent(m::MyDPM)\nupdate_suffstats!(m::MyDPM, data)\nupdate_suffstats!(m::MyDPM, data, i::Int, k0::Int, k1::Int)\nupdate_hyperpars!(rng::AbstractRNG, m::MyDPM, data)\nlogpredlik(m::MyDPM, data, i::Int, k::Int)","category":"page"},{"location":"create/","page":"Creating new DPMs","title":"Creating new DPMs","text":"see the documentation of each method for more details. Once again, we note that data can have any type, as long as each method is correctly implemented.","category":"page"},{"location":"start/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"start/#Using-a-custom-DPM","page":"Getting Started","title":"Using a custom DPM","text":"","category":"section"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Let MyDPM <: AbstractDPM be a datatype that conforms with the interface defined in this module, and let m an object of type MyDPM. Then, m contains the current state of the chain associated with the Gibbs sampler described in Neal's algorithm 3. Then, we can access the contents of m using the following accessors:","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"julia> DPMNeal3.n_clusters(m)      # return the number of active clusters\njulia> DPMNeal3.active_clusters(m) # return the active clusters\njulia> DPMNeal3.cluster_labels(m)  # return the cluster labels\njulia> DPMNeal3.cluster_sizes(m)   # return the cluster sizes\njulia> DPMNeal3.dp_mass(m)         # return the DP mass parameter","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Now, let MyData <: Any be the datatype of the sample expected by MyDPM, and let data an object of type MyData. Then, we can perform one Gibbs update using update!():","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"julia> rng = Random.MersenneTwister(1) # or any AbstractRNG object\njulia> update!(rng, m, data)","category":"page"},{"location":"start/","page":"Getting Started","title":"Getting Started","text":"Note that data can have any type, provided that MyDPM conforms with the interface.","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Install with the Julia package manager Pkg:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"# Press ']' to enter the Pkg REPL mode.\npkg> add DPMNeal3","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"julia> using Pkg; \njulia> Pkg.add(\"DPMNeal3\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DPMNeal3","category":"page"},{"location":"#DPMNeal3.jl","page":"Home","title":"DPMNeal3.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DPMNeal3.jl provides an interface for Dirichlet process mixture (DPM) models in Julia using Neal's algorithm 3 [1]. Particularly, this package provides:","category":"page"},{"location":"","page":"Home","title":"Home","text":"A method for performing one iteration of Neal's algorithm 3.\nSeveral methods for accessing the current state of the chain (e.g. the vector of cluster labels).","category":"page"}]
}
