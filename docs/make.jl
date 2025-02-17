import Pkg
Pkg.add("Documenter")
using Documenter
push!(LOAD_PATH, "../")
push!(LOAD_PATH, "../src/")
using MolSimToolkitShared
makedocs(
    modules=[MolSimToolkitShared],
    sitename="MolSimToolkitShared.jl",
    doctest=false,
    pages=[
        "Home" => "index.md",
    ],
)
deploydocs(
    repo="github.com/m3g/MolSimToolkitShared.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
