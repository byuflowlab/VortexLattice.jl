using Documenter, VortexLattice

makedocs(;
    format = Documenter.HTML(;size_threshold = 102400, size_threshold_warn = 102000),
    modules = [VortexLattice],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Advanced Usage" => "advanced.md",
        "Library" => "library.md",
        "Theory" => "theory.md"
    ],
    sitename = "VortexLattice.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu> and Andrew Ning <aning@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/VortexLattice.jl.git",
)
