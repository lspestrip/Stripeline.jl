if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, Stripeline

DocMeta.setdocmeta!(Stripeline, :DocTestSetup, :(using Stripeline); recursive=true)

ENV["GKSwstype"] = "100"  # Disable display of plots

makedocs(
    modules = [Stripeline],
    format = Documenter.HTML(
        assets = [
            Documenter.asset("https://threejs.org/build/three.js"),
            #joinpath("assets", "three.min.js"),
            joinpath("assets", "OrbitControls.js"),
        ],
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "Stripeline.jl",
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "Basic functions" => "basic.md",
        "Instrument database" => "instrumentdb.md",
        "Scanning strategy" => "scanning.md",
        "Data acquisition" => "acquisition.md",
        "Map making" => "mapmaking.md",
	"Full-scale simulation tutorial" => "simulation_tutorial.md",
    ],
    debug = true,
)

deploydocs(
    repo = "github.com/lspestrip/Stripeline.jl.git",
)
