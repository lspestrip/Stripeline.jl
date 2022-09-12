if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

library_path = joinpath("..", "src")
(library_path in LOAD_PATH) || push!(LOAD_PATH, library_path)

using Documenter, Stripeline

DocMeta.setdocmeta!(Stripeline, :DocTestSetup, :(using Stripeline); recursive=true)

ENV["GKSwstype"] = "100"  # Disable display of plots

makedocs(
    modules = [Stripeline],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "Stripeline.jl",
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "Basic functions" => "basic.md",
        "Instrument database" => "instrumentdb.md",
        "Time-ordered data" => "tod.md",
        "Scanning strategy" => "scanning.md",
        "Noise generation" => "noise.md",
        "Data acquisition" => "acquisition.md",
        "Map making" => "mapmaking.md",
	"Full-scale simulation tutorial" => "simulation_tutorial.md",
    ],
    debug = true,
)

deploydocs(
    repo = "github.com/lspestrip/Stripeline.jl.git",
)
