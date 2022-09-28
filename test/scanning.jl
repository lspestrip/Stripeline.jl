using Test
using Stripeline
using Dates
using AstroLib
using Healpix
#const Sl = Stripeline

# These are the old values that were used by @fincardona to test
# against the reference position for the Crab Nebula
const TEST_TENERIFE_LATITUDE_DEG = 28.3
const TEST_TENERIFE_LONGITUDE_DEG = -16.509722
const TEST_TENERIFE_HEIGHT_M = 2390

const eps = 3e-4 # Corresponds to pointing precision of 1.08 arcseconds

t_start = DateTime(2018, 01, 01, 0, 0, 0)

const crab_az_stellarium_rad = 3.168070267969909
const crab_alt_stellarium_rad = 1.4612458881566615

dirs = [π / 2 - crab_alt_stellarium_rad 2π - crab_az_stellarium_rad]

const crab_ra_astropy_rad = 1.4596726619436968
const crab_dec_astropy_rad = 0.3842255081802917
const crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

# Invert Crab coordinates into telescope pointing directions
groundq = telescopetoground(_ -> (0, deg2rad(20), 0), 0)
rotmatr = rotationmatrix_normalized(groundq)
vector = Float64[Healpix.ang2vec(dirs...)...]
dir = inv(rotmatr) * vector

# Compute skydirs
(skydirs, skyψ) = genpointings(_ -> (0, deg2rad(20), 0),
                               dir,
                               [0],
                               t_start,
                               latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
                               longitude_deg = TEST_TENERIFE_LONGITUDE_DEG,
                               height_m = TEST_TENERIFE_HEIGHT_M,
                               precession = true,
                               nutation = true,
                               aberration = true,
                               refraction = true)
crab_position_skydirs = sqrt((π/2 - skydirs[1])^2 + skydirs[2]^2)
@test skydirs[1] ≈ π/2 - crab_dec_astropy_rad atol = eps
@test skydirs[2] ≈ crab_ra_astropy_rad atol = eps
@test crab_position_skydirs ≈ crab_position atol = eps

# Check that Julian dates work too
(skydirs_jd, skyψ_jd) = genpointings(_ -> (0, deg2rad(20), 0),
                                     dir,
                                     [0],
                                     jdcnv(t_start),
                                     latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
                                     longitude_deg = TEST_TENERIFE_LONGITUDE_DEG,
                                     height_m = TEST_TENERIFE_HEIGHT_M,
                                     precession = true,
                                     nutation = true,
                                     aberration = true,
                                     refraction = true)

@test skydirs_jd[:] ≈ skydirs[:]
@test skyψ_jd[:] ≈ skyψ[:]

t_1 = DateTime(2019, 02, 01, 2, 0, 0)
t_2 = DateTime(2019, 02, 02, 2, 0, 0)
t_3 = DateTime(2019, 02, 03, 2, 0, 0)

days = [t_1, t_2, t_3]

dirs = [π / 2 - 0.6174703397894518 2π - 4.852881197778698 ;
        π / 2 - 0.6024799007695448 2π - 4.859519751514131 ;
        π / 2 - 0.5875049757874335 2π - 4.866139397516001]

# Invert crab coordinates into telescope pointing directions
skydirs = Array{Float64}(undef, 3, 2)
for (idx, day) in enumerate(days)
    local vector = Float64[Healpix.ang2vec(dirs[idx, 1], dirs[idx, 2])...]
    local dir = inv(rotmatr) * vector

    local (skydirections, skyψ) = genpointings(
        _ -> (0, deg2rad(20), 0),
        dir,
        [0],
        day,
        latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
        longitude_deg = TEST_TENERIFE_LONGITUDE_DEG,
        height_m = TEST_TENERIFE_HEIGHT_M,
        precession = true,
        nutation = true,
        aberration = true,
        refraction = true,
    )

    skydirs[idx, 1] = skydirections[1]
    skydirs[idx, 2] = skydirections[2]
end

crab_position_skydirs = sqrt.((π/2 .- skydirs[:, 1]).^2 + skydirs[:, 2].^2)

@test skydirs[1, 1] ≈ π/2 - crab_dec_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[2, 1] ≈ π/2 - crab_dec_astropy_rad atol = eps
@test skydirs[2, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[3, 1] ≈ π/2 - crab_dec_astropy_rad atol = eps
@test skydirs[3, 2] ≈ crab_ra_astropy_rad atol = eps
@test crab_position_skydirs[1] ≈ crab_position atol = eps
@test crab_position_skydirs[2] ≈ crab_position atol = eps
@test crab_position_skydirs[3] ≈ crab_position atol = eps

# Test flag `ground`
db = Sl.InstrumentDB()

time_duration = 86400
sampling_rate = 50.0
spin_velocity = 1

τ_s = 1 / sampling_rate
times = 0:τ_s:time_duration

(dirs, ψ) = genpointings(db.focalplane["I0"].orientation,
                         times;
                         latitude_deg = TEST_TENERIFE_LATITUDE_DEG) do time_s
    (0, deg2rad(20.0), Sl.timetorotang(time_s, spin_velocity))
end

expected_nsamples = convert(Int, time_duration * sampling_rate + 1)
@test size(dirs) == (expected_nsamples, 2)
@test size(ψ) == (expected_nsamples,)

(dirs, ψ) = genpointings(db.focalplane["I0"].orientation,
                         times;
                         ground = true,
                         latitude_deg = TEST_TENERIFE_LATITUDE_DEG) do time_s
    (0, deg2rad(20.0), Sl.timetorotang(time_s, spin_velocity))
end

expected_nsamples = convert(Int, time_duration * sampling_rate + 1)
@test size(dirs) == (expected_nsamples, 4)
@test size(ψ) == (expected_nsamples, 2)

# Check that all the values in the 3rd column are the same
@test dirs[1:end - 1, 3] ≈ dirs[2:end, 3]

# We strive for an accuracy of a few arcseconds
@test rad2deg(dirs[1, 3]) ≈ 20.0 atol = 1e-3

################################################################################

# Check that the sign of the wheel angles is correct

let defaultdb = InstrumentDB()
    wheelfn(time_s) = (0, deg2rad(20), Stripeline.timetorotang(time_s, 1))

    # We do not start from t = 0s because we want the longitude to be positive
    timerange = 7200:1:(7200 + 60)
    G0_vec = defaultdb.focalplane["G0"].orientation
    V0_vec = defaultdb.focalplane["V0"].orientation

    dirG0, psiG0 = Stripeline.genpointings(
        wheelfn,
        G0_vec,
        timerange,
        latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
    )
    dirV0, psiV0 = Stripeline.genpointings(
        wheelfn,
        V0_vec,
        timerange,
        latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
    )

    # We expect G0 to draw larger circles in the sky
    @test minimum(dirG0[:, 1]) < minimum(dirV0[:, 1])
    @test minimum(dirG0[:, 1]) < minimum(dirV0[:, 1])
    @test maximum(dirG0[:, 1]) > maximum(dirV0[:, 1])
    @test maximum(dirG0[:, 2]) > maximum(dirV0[:, 2])
end

#############################################################

# Test the PRM with non idealities
# Can't test wheel1ang_0 (the boresight motor zero point) because PyPRM doesn't support it, a solution must be found!

function angletomatrix(wheelanglesfn, time_s, config_ang::configuration_angles)
    rotationmatrix_normalized(telescopetoground(wheelanglesfn, time_s, config_ang))    
end

# Single configuration angles
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles()), 
                [0.9396926207859084 0.0 0.3420201433256687; 0.0 1.0 0.0; -0.3420201433256687 0.0 0.9396926207859084])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(forkang_rad=deg2rad(10))),
                [0.9396926207859084 0.0 0.3420201433256687; 0.0593911746138847 0.984807753012208 -0.16317591116653482; -0.33682408883346515 0.17364817766693033 0.9254165783983234])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(wheel2ang_0_rad=deg2rad(10))), 
                [0.984807753012208 0.0 0.17364817766693033; 0.0 1.0 0.0; -0.17364817766693033 0.0 0.984807753012208])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(wheel3ang_0_rad=deg2rad(-10))), 
                [0.9254165783983234 0.17364817766693033 0.33682408883346515; -0.16317591116653482 0.984807753012208 -0.0593911746138847; -0.3420201433256687 0.0 0.9396926207859084])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(zVAXang_rad=deg2rad(10))), 
                [0.9396926207859084 0.0 0.3420201433256687; 0.0593911746138847 0.984807753012208 -0.16317591116653482; -0.33682408883346515 0.17364817766693033 0.9254165783983234])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(omegaVAXang_rad=deg2rad(10))), 
                [0.9254165783983234 -0.17364817766693033 0.33682408883346515; 0.16317591116653482 0.984807753012208 0.0593911746138847; -0.3420201433256687 0.0 0.9396926207859084])

# Combination of different configuration angles
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles(wheel2ang_0_rad=deg2rad(48), wheel3ang_0_rad=deg2rad(-30))), 
                [0.7646550456261504 0.49999999999999994 -0.4065742997269626; -0.44147379642946344 0.8660254037844387 0.23473578139294538; 0.4694715627858908 0.0 0.882947592858927])
@test isapprox(angletomatrix(_ -> (0.0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles(wheel2ang_0_rad=deg2rad(48), wheel3ang_0_rad=deg2rad(-30), forkang_rad=deg2rad(52))),
				[0.882947592858927 2.95973511123774e-17 -0.4694715627858908; -0.36994863998783545 0.6156614753256583 -0.6957721980440043; 0.28903555496820393 0.788010753606722 0.5435968176547656])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles(wheel2ang_0_rad=deg2rad(48.),wheel3ang_0_rad=deg2rad(-30.),forkang_rad=deg2rad(52.),zVAXang_rad=deg2rad(42.))),
				[0.882947592858927 2.95973511123774e-17 -0.4694715627858908; -0.46832795365450275 -0.06975647374412532 -0.8807967768995136; -0.03274868074308757 0.9975640502598243 -0.06159131057870244])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles(wheel2ang_0_rad=deg2rad(48.),wheel3ang_0_rad=deg2rad(-30.),forkang_rad=deg2rad(52.),zVAXang_rad=deg2rad(42.),omegaVAXang_rad=deg2rad(73.))),
				[0.7060131423352384 0.06670844760071766 0.7050499456553592; 0.7074411401378279 -0.020394819144016727 -0.706477943457949; -0.03274868074308757 0.9975640502598243 -0.06159131057870244])
