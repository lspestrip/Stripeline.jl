using Test
using Stripeline

@testset "Instrument database" begin
    include("instrumentdb_tests.jl")
end
    
@testset "Scanning strategy" begin
    include("scanning.jl")
end

@testset "Map maker" begin
    include("map_maker_tests.jl")
end

@testset "TOD splitter" begin
    include("tod_splitter_tests.jl")
end

@testset "Noise generation" begin
    include("noisegeneration_tests.jl")
end
