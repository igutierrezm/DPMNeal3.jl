var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = DPMNeal3","category":"page"},{"location":"#DPMNeal3","page":"Home","title":"DPMNeal3","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DPMNeal3.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [DPMNeal3]","category":"page"},{"location":"#DPMNeal3.GenericBlock","page":"Home","title":"DPMNeal3.GenericBlock","text":"GenericBlock(rng::AbstractRNG, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)\n\nInitialize the generic block of a DPM, where N is the sample size, K0  is the initial number of clusters, and (a0, b0) determine the prior  distribution of the DP mass parameter, a Gamma(a0, 1 / b0)  distribution. \n\n\n\n\n\n","category":"type"},{"location":"#DPMNeal3.A-Tuple{GenericBlock}","page":"Home","title":"DPMNeal3.A","text":"Return the set of active clusters\n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.K-Tuple{GenericBlock}","page":"Home","title":"DPMNeal3.K","text":"Return the current number of cluster\n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.N-Tuple{GenericBlock}","page":"Home","title":"DPMNeal3.N","text":"Return the sample size\n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.P-Tuple{GenericBlock}","page":"Home","title":"DPMNeal3.P","text":"Return the set of passive clusters\n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.logpredlik-Tuple{Any, GenericBlock, Any, Any, Any}","page":"Home","title":"DPMNeal3.logpredlik","text":"logpredlik(sb, gb::GenericBlock, data, i, k)\n\nReturn the log(y_i  y_-i d_i d_-i). \n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.α-Tuple{GenericBlock}","page":"Home","title":"DPMNeal3.α","text":"Return the current number of cluster\n\n\n\n\n\n","category":"method"}]
}