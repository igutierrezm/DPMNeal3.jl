push!(LOAD_PATH,"../src/")
using DPMNeal3
using Documenter

DocMeta.setdocmeta!(DPMNeal3, :DocTestSetup, :(using DPMNeal3); recursive=true)

makedocs(;
    modules=[DPMNeal3],
    authors="Iván Gutiérrez <ivangutierrez1988@gmail.com> and contributors",
    repo="https://github.com/igutierrezm/DPMNeal3.jl/blob/{commit}{path}#{line}",
    sitename="DPMNeal3.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://igutierrezm.github.io/DPMNeal3.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "install.md",
        "The AbstractDPM type" => "abstractdpm.md",
        "The Skeleton type" => "skeleton.md",
        "Interface" => "interface.md",
    ],
)

deploydocs(;
    repo="github.com/igutierrezm/DPMNeal3.jl",
    devbranch="main",
)
