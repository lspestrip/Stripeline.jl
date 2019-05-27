if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, Stripeline

makedocs(
    modules = [Stripeline],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    sitename = "Stripeline.jl",
    pages = [
        "Introduction" => "index.md",
        "Basic functions" => "basic.md",
        "Instrument database" => "instrumentdb.md",
        "Scanning strategy" => "scanning.md",
        "Map making" => "mapmaking.md",
	"Simulation tutorial" => "simulation_tutorial.md",
    ])

deploydocs(
    repo = "github.com/lspestrip/Stripeline.jl.git",
)
