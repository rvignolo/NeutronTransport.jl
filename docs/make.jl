using NeutronTransport
using Documenter

makedocs(;
    modules=[NeutronTransport],
    authors="Ramiro Vignolo <ramirovignolo@gmail.com>",
    repo="https://github.com/rvignolo/NeutronTransport.jl/blob/{commit}{path}#L{line}",
    sitename="NeutronTransport.jl",
    format=Documenter.HTML(;
        prettyurls=false,
        canonical="https://rvignolo.github.io/FirstPkg.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rvignolo/NeutronTransport.jl",
)