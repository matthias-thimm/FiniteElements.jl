using FiniteElements
using Documenter

DocMeta.setdocmeta!(FiniteElements, :DocTestSetup, :(using FiniteElements); recursive=true)

makedocs(;
    modules=[FiniteElements],
    authors="Kai Partmann",
    sitename="FiniteElements.jl",
    format=Documenter.HTML(;
        canonical="https://kaipartmann.github.io/FiniteElements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kaipartmann/FiniteElements.jl",
    devbranch="main",
)
