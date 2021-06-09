var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = DPMNeal3","category":"page"},{"location":"#DPMNeal3","page":"Home","title":"DPMNeal3","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DPMNeal3.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [DPMNeal3]","category":"page"},{"location":"#DPMNeal3.GenericBlock","page":"Home","title":"DPMNeal3.GenericBlock","text":"GenericBlock(rng::AbstractRNG, N::Int; K0::Int = 1, a0 = 2.0, b0 = 4.0)\n\nInitialize the generic block of a DPM, where N is the sample size, K0 is  the initial number of clusters, and (a0, b0) determine the prior  distribution of the DP mass parameter, a Gamma(a0, 1 / b0)  distribution. \n\n\n\n\n\n","category":"type"},{"location":"#DPMNeal3.logh-Tuple{Any, GenericBlock, Any, Any}","page":"Home","title":"DPMNeal3.logh","text":"logh(sb, gb::GenericBlock, data, i)\n\nReturn the log of h(y_i) = int q(y_i  theta) g_0(theta) dtheta\n\n\n\n\n\n","category":"method"},{"location":"#DPMNeal3.logq-Tuple{Any, GenericBlock, Any, Any, Any}","page":"Home","title":"DPMNeal3.logq","text":"logq(sb, gb::GenericBlock, data, i)\n\nReturn the log of q(y_i) = p(y_i  d_i = k d_-i)\n\n\n\n\n\n","category":"method"}]
}
