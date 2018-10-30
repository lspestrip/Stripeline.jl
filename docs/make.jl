push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, Stripeline

makedocs(modules = [Stripeline],
    format = :html,
    sitename = "Stripeline.jl",
    pages = Any[
        "Introduction" => "index.md",
        "Instrument database" => "instrumentdb.md",
        "Scanning strategy" => "scanning.md",
	"Simulation tutorial" => "simulation_tutorial.md"
    ])

deploydocs(repo = "github.com/lspestrip/Stripeline.jl.git",
    target = "build",
    julia = "0.7",
    deps = nothing,
    make = nothing)
