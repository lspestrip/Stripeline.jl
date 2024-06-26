using Test
using Stripeline
using Dates
using AstroLib
using Healpix
using StaticArrays
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
(skydirs, skyψ) = genpointings(
    _ -> (0, deg2rad(20), 0),
    dir,
    [0],
    t_start,
    latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
    longitude_deg = TEST_TENERIFE_LONGITUDE_DEG,
    height_m = TEST_TENERIFE_HEIGHT_M,
    precession = true,
    nutation = true,
    aberration = true,
    refraction = true,
)
crab_position_skydirs = sqrt((π / 2 - skydirs[1])^2 + skydirs[2]^2)
@test skydirs[1] ≈ π / 2 - crab_dec_astropy_rad atol = eps
@test skydirs[2] ≈ crab_ra_astropy_rad atol = eps
@test crab_position_skydirs ≈ crab_position atol = eps

# Check that Julian dates work too
(skydirs_jd, skyψ_jd) = genpointings(
    _ -> (0, deg2rad(20), 0),
    dir,
    [0],
    jdcnv(t_start),
    latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
    longitude_deg = TEST_TENERIFE_LONGITUDE_DEG,
    height_m = TEST_TENERIFE_HEIGHT_M,
    precession = true,
    nutation = true,
    aberration = true,
    refraction = true,
)

@test skydirs_jd[:] ≈ skydirs[:]
@test skyψ_jd[:] ≈ skyψ[:]

t_1 = DateTime(2019, 02, 01, 2, 0, 0)
t_2 = DateTime(2019, 02, 02, 2, 0, 0)
t_3 = DateTime(2019, 02, 03, 2, 0, 0)

days = [t_1, t_2, t_3]

dirs = [
    π/2-0.6174703397894518 2π-4.852881197778698
    π/2-0.6024799007695448 2π-4.859519751514131
    π/2-0.5875049757874335 2π-4.866139397516001
]

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

crab_position_skydirs = sqrt.((π / 2 .- skydirs[:, 1]) .^ 2 + skydirs[:, 2] .^ 2)

@test skydirs[1, 1] ≈ π / 2 - crab_dec_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[2, 1] ≈ π / 2 - crab_dec_astropy_rad atol = eps
@test skydirs[2, 2] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[3, 1] ≈ π / 2 - crab_dec_astropy_rad atol = eps
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

(dirs, ψ) = genpointings(
    db.focalplane["I0"].orientation,
    times;
    latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
) do time_s
    (0, deg2rad(20.0), Sl.timetorotang(time_s, spin_velocity))
end

expected_nsamples = convert(Int, time_duration * sampling_rate + 1)
@test size(dirs) == (expected_nsamples, 2)
@test size(ψ) == (expected_nsamples,)

(dirs, ψ) = genpointings(
    db.focalplane["I0"].orientation,
    times;
    ground = true,
    latitude_deg = TEST_TENERIFE_LATITUDE_DEG,
) do time_s
    (0, deg2rad(20.0), Sl.timetorotang(time_s, spin_velocity))
end

expected_nsamples = convert(Int, time_duration * sampling_rate + 1)
@test size(dirs) == (expected_nsamples, 4)
@test size(ψ) == (expected_nsamples, 2)

# Check that all the values in the 3rd column are the same
@test dirs[1:end-1, 3] ≈ dirs[2:end, 3]

# We strive for an accuracy of a few arcseconds
@test rad2deg(dirs[1, 3]) ≈ 20.0 atol = 1e-3

################################################################################

# Check that the sign of the wheel angles is correct

let defaultdb = InstrumentDB()
    wheelfn(time_s) = (0, deg2rad(20), Stripeline.timetorotang(time_s, 1))

    # We do not start from t = 0s because we want the longitude to be positive
    timerange = 7200:1:(7200+60)
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
@test isapprox(
    angletomatrix(_ -> (0.0, deg2rad(20.0), 0), 0, configuration_angles()),
    [
        0.9396926207859084 0.0 0.3420201433256687
        0.0 1.0 0.0
        -0.3420201433256687 0.0 0.9396926207859084
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(forkang_rad = deg2rad(10)),
    ),
    [
        0.9396926207859084 0.0 0.3420201433256687
        0.0593911746138847 0.984807753012208 -0.16317591116653482
        -0.33682408883346515 0.17364817766693033 0.9254165783983234
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(wheel2ang_0_rad = deg2rad(10)),
    ),
    [
        0.984807753012208 0.0 0.17364817766693033
        0.0 1.0 0.0
        -0.17364817766693033 0.0 0.984807753012208
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(wheel3ang_0_rad = deg2rad(-10)),
    ),
    [
        0.9254165783983234 0.17364817766693033 0.33682408883346515
        -0.16317591116653482 0.984807753012208 -0.0593911746138847
        -0.3420201433256687 0.0 0.9396926207859084
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(zVAXang_rad = deg2rad(10)),
    ),
    [
        0.9396926207859084 0.0 0.3420201433256687
        0.0593911746138847 0.984807753012208 -0.16317591116653482
        -0.33682408883346515 0.17364817766693033 0.9254165783983234
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(omegaVAXang_rad = deg2rad(10)),
    ),
    [
        0.9254165783983234 -0.17364817766693033 0.33682408883346515
        0.16317591116653482 0.984807753012208 0.0593911746138847
        -0.3420201433256687 0.0 0.9396926207859084
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(rollang_rad = deg2rad(45)),
    ),
    [
        0.2218884684027577 -0.928995249589305 0.29619813272602386
        0.9446039478901318 0.28014092350145736 0.17101007166283433
        -0.24184476264797528 0.24184476264797522 0.9396926207859084
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(panang_rad = deg2rad(25)),
    ),
    [
        0.8137976813493738 -0.32797515353481177 0.4797558050656603
        0.46984631039295416 0.8571575464476954 -0.21101039116094453
        -0.3420201433256687 0.39713126196710286 0.8516507396391465
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(tiltang_rad = deg2rad(10)),
    ),
    [
        0.7500000000000001 -0.49999999999999994 0.4330127018922193
        0.4330127018922193 0.8660254037844387 0.24999999999999994
        -0.49999999999999994 0.0 0.8660254037844387
    ],
)

# Combination of different configuration angles
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), 0),
        0,
        configuration_angles(wheel2ang_0_rad = deg2rad(48), wheel3ang_0_rad = deg2rad(-30)),
    ),
    [
        0.7646550456261504 0.49999999999999994 -0.4065742997269626
        -0.44147379642946344 0.8660254037844387 0.23473578139294538
        0.4694715627858908 0.0 0.882947592858927
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0.0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(
            wheel2ang_0_rad = deg2rad(48),
            wheel3ang_0_rad = deg2rad(-30),
            forkang_rad = deg2rad(52),
        ),
    ),
    [
        0.882947592858927 2.95973511123774e-17 -0.4694715627858908
        -0.36994863998783545 0.6156614753256583 -0.6957721980440043
        0.28903555496820393 0.788010753606722 0.5435968176547656
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(
            wheel2ang_0_rad = deg2rad(48.0),
            wheel3ang_0_rad = deg2rad(-30.0),
            forkang_rad = deg2rad(52.0),
            zVAXang_rad = deg2rad(42.0),
        ),
    ),
    [
        0.882947592858927 2.95973511123774e-17 -0.4694715627858908
        -0.46832795365450275 -0.06975647374412532 -0.8807967768995136
        -0.03274868074308757 0.9975640502598243 -0.06159131057870244
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(
            wheel2ang_0_rad = deg2rad(48.0),
            wheel3ang_0_rad = deg2rad(-30.0),
            forkang_rad = deg2rad(52.0),
            zVAXang_rad = deg2rad(42.0),
            omegaVAXang_rad = deg2rad(73.0),
        ),
    ),
    [
        0.7060131423352384 0.06670844760071766 0.7050499456553592
        0.7074411401378279 -0.020394819144016727 -0.706477943457949
        -0.03274868074308757 0.9975640502598243 -0.06159131057870244
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(
            rollang_rad = deg2rad(68),
            tiltang_rad = deg2rad(91),
            panang_rad = deg2rad(137),
        ),
    ),
    [
        0.47444246561235864 0.34112613016341686 0.8115031177656669
        -0.21412297968452368 -0.8494536195128901 0.4822653811621473
        0.8538475838196693 -0.4025686421173186 -0.3299739262261965
    ],
)
@test isapprox(
    angletomatrix(
        _ -> (0, deg2rad(20.0), deg2rad(-30.0)),
        0,
        configuration_angles(
            wheel2ang_0_rad = deg2rad(48.0),
            wheel3ang_0_rad = deg2rad(-30.0),
            forkang_rad = deg2rad(52.0),
            zVAXang_rad = deg2rad(42.0),
            omegaVAXang_rad = deg2rad(73.0),
            rollang_rad = deg2rad(26),
            panang_rad = deg2rad(33),
            tiltang_rad = deg2rad(79),
        ),
    ),
    [
        -0.17570296000187266 0.5751787563280825 0.7989354592928397
        0.45810521448987807 -0.670565309868962 0.5835081641738409
        0.8713599040027968 0.4685206115735181 -0.14567207771914797
    ],
)

# New test to check that the polarization angle is computed correctly
# if Julian Dates are used

let times = 0.0:0.1:1.0, latitude_deg = 28.29
    (_, psi_solar) = Sl.genpointings(
        SVector(0.0, 0.0, 1.0),
        times,
        Dates.DateTime(2026, 01, 01, 00, 00, 00),
        latitude_deg=latitude_deg,
    ) do time_s
        (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1.0))
    end

    @test psi_solar[1] ≈ 1.5747867612568307
    @test psi_solar[2] ≈ 1.560929162843191
    @test psi_solar[3] ≈ 1.5470723919649034
    @test psi_solar[4] ≈ 1.533217620642667
    @test psi_solar[5] ≈ 1.5193660189782203
    @test psi_solar[6] ≈ 1.5055187540549246
    @test psi_solar[7] ≈ 1.4916769888231458
    @test psi_solar[8] ≈ 1.4778418810155476
    @test psi_solar[9] ≈ 1.4640145820540134
    @test psi_solar[10] ≈ 1.4501962359892084
    @test psi_solar[11] ≈ 1.4363879784382272
end
