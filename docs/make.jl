using Documenter, VortexFlow

makedocs(;
    modules = [VortexFlow],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Library" => "library.md",
        "Theory" => "theory.md"
    ],
    sitename = "VortexFlow.jl",
    authors = ["Taylor McDonnell <taylormcd@byu.edu>", "Andrew Ning <aning@byu.edu>"],
)

deploydocs(
    repo = "github.com/byuflowlab/VortexFlow.jl.git",
)
