using Documenter
using Pkg
Pkg.activate(".") 
Pkg.instantiate()
using ATA

makedocs(
    sitename="ATA",
    format=Documenter.HTML(),
    modules=[ATA],
    doctest=true,
    pages = [
        "index.md",
        "utils.md",
        "build.md",
        "opt.md",
        "examples.md" => [
            "jump_example.md",
            "siman_example.md",
            "custom_example.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/giadasp/ATA.jl.git",
    devurl = "docs"
)
