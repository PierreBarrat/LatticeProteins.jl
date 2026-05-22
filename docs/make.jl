using LatticeProteins
using Documenter

DocMeta.setdocmeta!(LatticeProteins, :DocTestSetup, :(using LatticeProteins); recursive=true)

makedocs(;
    modules=[LatticeProteins],
    authors="Pierre Barrat-Charlaix",
    sitename="LatticeProteins.jl",
    format=Documenter.HTML(;
        canonical="https://pierrebarrat.github.io/LatticeProteins.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Index" => "index.md", "Reference" => "reference.md"],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/PierreBarrat/LatticeProteins.jl.git", devbranch="main")
