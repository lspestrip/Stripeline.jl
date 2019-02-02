using Quaternions
import Healpix
using StaticArrays
using LinearAlgebra
using AstroLib
using Dates

export TENERIFE_LATITUDE_DEG, TENERIFE_LONGITUDE_DEG, TENERIFE_HEIGHT_M
export timetorotang, telescopetoground, groundtoearth
export genpointings, polarizationangle

const TENERIFE_LATITUDE_DEG = 28.3
const TENERIFE_LONGITUDE_DEG = -16.509722
const TENERIFE_HEIGHT_M = 2390


"""
    timetorotang(time, rpm)

Convert a time into a rotation angle, given the number of rotations per minute.
The time should be expressed in seconds. The return value is in radians.
`time` can either be a scalar or a vector.
"""
function timetorotang(time_s, rpm)
    if rpm == 0
        0.0
    else
        2 * π * time_s * (rpm / 60)
    end
end


"""
    telescopetoground(wheelanglesfn, time_s)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transform from the focal plane to the ground of the
telescope. The parameter `wheelanglesfn` must be a function which
takes as input a time, `time_s`, in seconds, and it must return a
3-tuple containing the angles of the following motors:

1. The boresight motor
2. The altitude motor
3. The ground motor

Example:
`````julia
telescopetoground(3600.0) do
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````
"""
function telescopetoground(wheelanglesfn, time_s)
    (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)
    
    qwheel1 = qrotation([0, 0, 1], wheel1ang)
    qwheel2 = qrotation([1, 0, 0], wheel2ang)
    qwheel3 = qrotation([0, 0, 1], wheel3ang)
    
    qwheel3 * (qwheel2 * qwheel1)
end


"""
    groundtoearth(groundq, time_s, latitude_deg; day_duration_s=86400.0)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transformation from the ground of the telescope to the
Equatorial coordinate system. The parameter `groundq` must be a
quaternion describing the coordinate transformation from the focal
plane of the telescope to the ground. The parameter `time_s` must be a
time in seconds, and `latitude_deg` is the latitude (in degrees, N is
positive) of the location where the observation is made.

The keyword `day_duration_s` specifies the length of a day in seconds.

"""
function groundtoearth(groundq, time_s, latitude_deg; day_duration_s=86400.0)
    locq = qrotation([1, 0, 0], deg2rad(90 - latitude_deg))
    earthq = qrotation([0, 0, 1], 2 * π * time_s / day_duration_s)
    
    earthq * (locq * groundq)
end


"""
    vector2equatorial(dir, jd, latitude_deg, longitude_deg, height_m)

Transform the Healpix coordinates of a vector into equatorial coordinates. 
The parameter `vector` is 3D vector, `jd` is the julian date. The 
paramters `latitude_deg`, `longitude_deg` and `height_m` should contain the 
latitude (in degrees, N is positive), the longitude (in degrees, counterclockwise
is positive) and the height (in meters) of the location where the observation is 
made. 
"""
function vector2equatorial(vector, jd, latitude_deg, longitude_deg, height_m)
    (θ, ϕ) = Healpix.vec2ang(vector[1], vector[2], vector[3])
    Alt_rad = π/2 - θ 
    Az_rad = 2π - ϕ
    
    Ra_deg, Dec_deg, HA_deg = AstroLib.hor2eq(rad2deg(Alt_rad),
                                              rad2deg(Az_rad),
                                              jd,
                                              latitude_deg,
                                              longitude_deg,
                                              height_m,
                                              precession=true,
                                              nutate=true,
                                              aberration=true)
    (deg2rad(Dec_deg), deg2rad(Ra_deg))
end


"""
    polarizationangle(northdir, poldir)

Calculate the polarization angle projected in the sky in IAU conventions.
The parameter `northdir` must be a versor that points the North and `poldir`
must be a versor that identify the polarization direction projected in the sky.
The return value is in radians.
"""
function polarizationangle(northdir, poldir)
    cosψ = clamp(dot(northdir, poldir), -1, 1)
    crosspr = northdir × poldir
    sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
    ψ = atan(sinψ, cosψ) 
    ψ
end


"""
    genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0, 
                 ground=false)

Generate a set of pointings for some STRIP detector. The parameter
`wheelanglesfn` must be a function which takes as input a time in seconds
and returns a 3-tuple containing the angles (in radians) of the three
motors:
1. The boresight motor
2. The altitude motor
3. The ground motor

The parameter `dir` must be a normalized vector which tells the pointing
direction of the beam (boresight is [0, 0, 1]). The parameter `timerange_s`
is either a range or a vector which specifies at which times (in second)
the pointings should be computed. The keyword `latitude_deg` should contain
the latitude (in degrees, N is positive) of the location where the observation
is made. The keyword `ground` must be a boolean: if true the angles will be 
referred to the ground coordinate system otherwise they will be expressed in 
equatorial coordinates; default is false.

Return a 2-tuple containing the directions (a N×2 array containing the
colatitude and the longitude) and the polarization angles at each time step.

Example:
`````julia
genpointings([0, 0, 1], 0:0.1:1) do time_s
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    return (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````
"""
function genpointings(wheelanglesfn,
                      dir,
                      timerange_s;
                      latitude_deg=0.0,
                      ground=false)
    
    dirs = Array{Float64}(undef, length(timerange_s), 2)
    ψ = Array{Float64}(undef, length(timerange_s))

    zaxis = [1; 0; 0] #should be different for each horn...
    for (idx, time_s) = enumerate(timerange_s)
        
        # This is in the ground reference frame
        groundq = telescopetoground(wheelanglesfn, time_s)
        
        if ground
            rotmatr = rotationmatrix(groundq)
        else
            # Now from the ground reference frame to the Earth reference frame
            quat = groundtoearth(groundq, time_s, latitude_deg)

            rotmatr = rotationmatrix(quat)
        end
        
        vector = rotmatr * dir
        poldir = rotmatr * zaxis

        # The North for a vector v is just -dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(vector[1], vector[2], vector[3])
        dirs[idx, 1] = θ
        dirs[idx, 2] = ϕ
        northdir = @SArray [-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ)]
        ψ[idx] = polarizationangle(northdir, poldir)
    end
    
    (dirs, ψ)
end


"""
    genpointings(wheelanglesfn, dir, timerange_s, t_start, t_stop; 
                 latitude_deg=0.0, longitude_deg=0.0, height_m=0.0)

Generate a set of pointings for some STRIP detector. The parameter
`wheelanglesfn` must be a function which takes as input a time in seconds
and returns a 3-tuple containing the angles (in radians) of the three
motors:
1. The boresight motor
2. The altitude motor
3. The ground motor

The parameter `dir` must be a normalized vector which tells the pointing
direction of the beam (boresight is [0, 0, 1]). The parameter `timerange_s`
is either a range or a vector which specifies at which times (in second)
the pointings should be computed. The parameter `t_start` and `t_start` must be 
two DateTime which tell the exact UTC date and time of the observation. The 
keywords `latitude_deg`, `longitude_deg` and `height_m` should contain the 
latitude (in degrees, N is positive), the longitude (in degrees, counterclockwise
is positive) and the height (in meters) of the location where the observation is 
made.

Return a 2-tuple containing the sky directions (a N×2 array containing the 
Declination and the RightAscension) and the polarization angle given in 
equatorial coordinates at each time step.

Example:
`````julia
using Dates

genpointings([0, 0, 1], 
             0:0.1:1, 
             DateTime(2019, 01, 01, 0, 0, 0), 
             DateTime(2022, 04, 13, 21, 10, 10), 
             latitude_deg=10.0,
             longitude_deg=20.0,
             height_m=1000) do time_s
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    return (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````
"""
function genpointings(wheelanglesfn,
                      dir,
                      timerange_s,
                      t_start,
                      t_stop;
                      latitude_deg=0.0,
                      longitude_deg=0.0,
                      height_m=0)
    
    skydirs = Array{Float64}(undef, length(timerange_s), 2)
    skyψ = Array{Float64}(undef, length(timerange_s))
    
    jd_start = AstroLib.jdcnv(t_start)
    jd_stop = AstroLib.jdcnv(t_stop)
    jd_range = range(jd_start, stop=jd_stop, length=length(timerange_s))

#    zaxis = [1; 0; 0] #should be different for each horn...
    for (idx, time_s) = enumerate(timerange_s)

        # This is in the ground reference frame
        groundq = telescopetoground(wheelanglesfn, time_s)
        
        rotmatr = rotationmatrix(groundq)
        vector = rotmatr * dir

        Dec_rad, Ra_rad = vector2equatorial(vector,
                                            jd_range[idx],
                                            latitude_deg,
                                            longitude_deg,
                                            height_m)

        skydirs[idx, 1] = Dec_rad
        skydirs[idx, 2] = Ra_rad

        # poldir = rotmatr * zaxis
        # Decpol_rad, Rapol_rad = vector2equatorial(poldir,
        #                                           jd_range[idx],
        #                                           latitude_deg,
        #                                           longitude_deg,
        #                                           height_m)
        
        # skypoldir = @SArray [cos(Decpol_rad) * cos(Rapol_rad),
        #                      cos(Decpol_rad) * sin(Rapol_rad),
        #                      sin(Decpol_rad)]
        # skynorthdir = @SArray [-sin(Decpol_rad) * cos(Rapol_rad),
        #                        -sin(Decpol_rad) * sin(Rapol_rad),
        #                        cos(Decpol_rad)]
        
        # skyψ[idx] = polarizationangle(skynorthdir, skypoldir)
        skyψ[idx] = 0.
    end

    (skydirs, skyψ) # The polarization angle is still missing
end
