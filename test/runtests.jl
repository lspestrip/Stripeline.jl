using Test
using Stripeline
using Documenter

@testset "Doctests" begin
    DocMeta.setdocmeta!(Stripeline, :DocTestSetup, :(using Stripeline); recursive=true)
    doctest(Stripeline, manual = false)
end

@testset "Instrument database" begin
    include("instrumentdb_tests.jl")
end

@testset "Time-oriented data" begin
    include("tod_tests.jl")
end

@testset "Scanning strategy" begin
    include("scanning.jl")
end

@testset "Pointing reconstruction star tracker" begin
    include("scanning_st.jl")
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

@testset "ADC simulation" begin
    include("adc_tests.jl")
end
