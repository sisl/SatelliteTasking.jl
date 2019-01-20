using Documenter, SatelliteTasking

include("src/makeplots.jl")

makedocs(
    modules   = [SatelliteTasking],
    doctest   = false,
    clean     = true,
    linkcheck = true,
    format    = Documenter.HTML(),
    sitename  = "SatelliteTasking.jl",
    authors   = "Duncan Eddy",
    pages     = Any[
        "Home" => "index.md",
        "Modules" => Any[
            "Data Structures" => "modules/data_structures.md",
        ],
        "Function Index" => "function_index.md",
    ]
)

deploydocs(
    repo = "github.com/sisl/SatelliteTasking.jl",
    devbranch = "master",
    devurl = "latest",
    deps = makeplots,
)