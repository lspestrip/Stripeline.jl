push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, Stripeline2

makedocs(modules = [Stripeline2],
    format = :html,
    sitename = "Stripeline2.jl",
    pages = Any[
        "Introduction" => "index.md",
        "Instrument database" => "instrumentdb.md",
        "Random numbers" => "rng.md",
        "Scanning strategy" => "scanning.md",
    ])

deploydocs(repo = "github.com/ziotom78/Stripeline2.jl.git",
    target = "build",
    julia = "0.6",
    deps = nothing,
    make = nothing)
