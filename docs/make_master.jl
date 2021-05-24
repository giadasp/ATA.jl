using Documenter
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using ATA

makedocs(
    sitename = "ATA",
    format = Documenter.HTML(),
    modules = [ATA],
    doctest = true,
    pages = ["index.md", "model_build.md", "opt.md", "compact.md", "print_and_plot.md", "structs.md", "utils.md", "examples.md"],
)

deploydocs(repo = "github.com/giadasp/ATA.jl.git", devurl = "docs", devbranch = "master")
