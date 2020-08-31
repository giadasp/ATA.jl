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
        "opt.md"
    ]
)

deploydocs(
    repo = "github.com/giadasp/ATA.jl.git"
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#= deploydocs(
    repo = "<repository url>"
) =#
