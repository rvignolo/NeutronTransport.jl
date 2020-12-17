using NeutronTransport
using Documenter

makedocs(;
    modules=[NeutronTransport],
    authors="Ramiro Vignolo <ramirovignolo@gmail.com>",
    repo="https://github.com/rvignolo/NeutronTransport.jl/blob/{commit}{path}#L{line}",
    sitename="NeutronTransport.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
