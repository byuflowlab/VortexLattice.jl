using Documenter, VLMFlow

makedocs(;
    modules = [VLMFlow],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Library" => "library.md",
        "Theory" => "theory.md"
    ],
    sitename = "VLMFlow.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu> and Andrew Ning <aning@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/VLMFlow.jl.git",
)
