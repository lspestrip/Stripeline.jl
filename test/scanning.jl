using Test
using Stripeline
using Dates
using AstroLib
using Healpix
const Sl = Stripeline

# These are the old values that were used by @fincardona to test
# against the reference position for the Crab Nebula
const TEST_TENERIFE_LATITUDE_DEG = 28.3
const TEST_TENERIFE_LONGITUDE_DEG = -16.509722
const TEST_TENERIFE_HEIGHT_M = 2390

eps = 3e-4 # Corresponds to pointing precision of 1.08 arcseconds

t_start = DateTime(2018, 01, 01, 0, 0, 0)

crab_az_stellarium_rad = 3.168070267969909
crab_alt_stellarium_rad = 1.4612458881566615

dirs = [π / 2 - crab_alt_stellarium_rad 2π - crab_az_stellarium_rad]

crab_ra_astropy_rad = 1.4596726619436968
crab_dec_astropy_rad = 0.3842255081802917
crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

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
crab_position_skydirs = sqrt(skydirs[1]^2 + skydirs[2]^2)
@test skydirs[1] ≈ crab_dec_astropy_rad atol = eps
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

crab_ra_astropy_rad = 1.4596726619436968
crab_dec_astropy_rad = 0.3842255081802917
crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

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

crab_position_skydirs = sqrt.(skydirs[:, 1].^2 + skydirs[:, 2].^2)

@test skydirs[1, 1] ≈ crab_dec_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[2, 1] ≈ crab_dec_astropy_rad atol = eps
@test skydirs[2, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[3, 1] ≈ crab_dec_astropy_rad atol = eps
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
    quaternion = telescopetoground(wheelanglesfn, time_s, config_ang)
    rotationmatrix_normalized(quaternion)    
end

# All the angles set to 0 and no "non idealities" simply test to see if angletomatrix and configuration_angles declaration work well
config = configuration_angles()
@test angletomatrix(_ -> (0, 0, 0), 0, config) == [1.0 0 0; 0 1.0 0; 0 0 1.0]

# Ideal
@test angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, config) == [0.9396926207859084 0.0 0.3420201433256687; 0.0 1.0 0.0; -0.3420201433256687 0.0 0.9396926207859084]

# tiltFork_10 deg
# config = configuration_angles(forkang = deg2rad(10))
# @test angletomatrix(_ -> (0.0, 20.0, 0), 0, config) == [0.9396926207859084 0.0593911746138847 -0.33682408883346515; 0.0 0.984807753012208 0.17364817766693033; 0.9396926207859084 0.984807753012208 0.9254165783983234]

# ctheta0_10 deg
# config = configuration_angles(wheel2ang_0 = deg2rad(10))
# @test angletomatrix(_ -> (0.0, 20.0, 0), 0, config) == [0.984807753012208 0.0 -0.17364817766693033; 0.0 1.0 0.0; 0.984807753012208 1.0 0.984807753012208]

# cphi0_10 deg
# config = configuration_angles(wheel3ang_0 = deg2rad(10))
# @test angletomatrix(_ -> (0.0, 20.0, 0), 0, config) == [0.9254165783983234 -0.16317591116653482 -0.3420201433256687; 0.17364817766693033 0.984807753012208 0.0; 0.9254165783983234 0.984807753012208 0.9396926207859084]

# zVAX_10 deg
# config = configuration_angles(zVAXang = deg2rad(10))
# @test angletomatrix(_ -> (0.0, 20.0, 0), 0, config) == [0.9396926207859084 0.0593911746138847 -0.33682408883346515; 0.0 0.984807753012208 0.17364817766693033; 0.9396926207859084 0.984807753012208 0.9254165783983234]

# omegaVAX_10 deg
# config = configuration_angles(omegaVAXang = deg2rad(10))
# @test angletomatrix(_ -> (0.0, 20.0, 0), 0, config) == [0.9254165783983234 0.16317591116653482 -0.3420201433256687; -0.17364817766693033 0.984807753012208 0.0; 0.9254165783983234 0.984807753012208 0.9396926207859084]