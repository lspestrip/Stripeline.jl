using Test
using Stripeline

# Scanning strategy

(dirs, psi) = genpointings([0, 0, 1], 0:20:60, latitude_deg=0.0) do time_s
    wheel1ang = 0.0
    wheel2ang = 30.0
    wheel3ang = timetorotang(time_s, 1.0)

    (wheel1ang, wheel2ang, wheel3ang)
end

@test size(dirs) == (4, 2)
@test length(psi) == 4
#@test dirs ≈ [0.15487 4.71239; 2.0875 3.3214; 2.0875 6.10774; 0.15487 4.71675], atol = 1e-5
# This fails, but we have to revise the way polarization angles are calculated
#@test psi ≈ [0, π / 2, π, 3π / 4]

# Instrument DB

include("instrumentdb_tests.jl")

# Map-maker

include("map_maker_tests.jl")

# TOD-splitter

include("tod_splitter_tests.jl")

# Scanning 

include("scanning.jl")
