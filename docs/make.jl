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
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "start.md",
        "Create new DPM models" => "create.md",
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/igutierrezm/DPMNeal3.jl",
)
