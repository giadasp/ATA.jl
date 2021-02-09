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
    pages = ["index.md", "utils.md", "build.md", "opt.md", "print.md", "examples.md"],
)

deploydocs(repo = "github.com/giadasp/ATA.jl.git", devurl = "docs", devbranch = "release")
