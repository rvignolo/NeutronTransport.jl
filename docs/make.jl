using NeutronTransport
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

makedocs(bib;
    modules=[NeutronTransport],
    authors="Ramiro Vignolo <ramirovignolo@gmail.com>",
    repo="https://github.com/rvignolo/NeutronTransport.jl/blob/{commit}{path}#L{line}",
    sitename="NeutronTransport.jl",
    format=Documenter.HTML(;
        prettyurls=false,
        canonical="https://rvignolo.github.io/NeutronTransport.jl",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Method of Characteristics" => "moc.md",
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/rvignolo/NeutronTransport.jl",
    devbranch = "main"
)