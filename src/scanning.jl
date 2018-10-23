using Quaternions
import Healpix
using StaticArrays
using LinearAlgebra

export TENERIFE_LATITUDE_DEG, TENERIFE_LONGITUDE_DEG, TENERIFE_HEIGHT_M
export timetorotang, genpointings, genskypointings, genastropyskypointings

TENERIFE_LATITUDE_DEG = 28.3
TENERIFE_LONGITUDE_DEG = -16.509722
TENERIFE_HEIGHT_M = 2390

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
    genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)

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
is made.

Return a 2-tuple containing the directions (a N×2 array containing the
colatitude and the longitude) and the polarization angles at each time step.
Directions are expressed in equatorial coordinates.

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
function genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)
    
    dirs = Array{Float64}(undef, length(timerange_s), 2)
    ψ = Array{Float64}(undef, length(timerange_s))

    zaxis = [1; 0; 0]
    for (idx, time_s) = enumerate(timerange_s)
        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)
        
        qwheel1 = qrotation([0, 0, 1], wheel1ang)
        qwheel2 = qrotation([1, 0, 0], wheel2ang)
        qwheel3 = qrotation([0, 0, 1], wheel3ang)
        
        # This is in the ground reference frame
        groundq = qwheel3 * (qwheel2 * qwheel1)
        
        # Now from the ground reference frame to the Earth reference frame
        locq = qrotation([1, 0, 0], deg2rad(90 - latitude_deg))
        earthq = qrotation([0, 0, 1], 2 * π * time_s / 86400)

        quat = earthq * (locq * groundq)
        rotmatr = rotationmatrix(quat)
        
        vector = rotmatr * dir
        poldir = rotmatr * zaxis

        # The North for a vector v is just -dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(vector[1], vector[2], vector[3])
        dirs[idx, 1] = θ
        dirs[idx, 2] = ϕ
        northdir = @SArray [-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ)]
        
        cosψ = clamp(dot(northdir, poldir), -1, 1)
        crosspr = northdir × poldir
        sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
        ψ[idx] = atan(cosψ, sinψ)
    end
    
    (dirs, ψ)
end


"""
    genastropyskypointings(t_start, t_stop, dirs; latitude_deg=0.0, 
                           longitude_deg=0.0, height_m=0.0)

Project a set of pointings in the sky exploiting `astropy`. Very slow but more 
accurate.  

The parameter `t_start` and `t_start` must be two strings which tells the exact 
UTC date and time of the observation in "iso" format. The parameter `dirs` must 
be a N×2 array containing the observed directions expressed in colatitude and 
the longitude. The keywords `latitude_deg`, `longitude_deg` and `height_m` should
contain the latitude (in degrees, N is positive), the longitude (in degrees, 
counterclockwise is positive) and the height (in meters)  of the location where 
the observation is made.

Return a 2-tuple containing the observed directions projected in sky coordinates 
(a N×2 array containing the Right ascension and the Declination, in radians) at 
each time step. 
Directions are expressed in ICRS coordinates.

Example:
`````julia
genskypointings("2019-01-01 00:00:00", "2020-01-25 14:52:10.05", dirs)
`````
"""
function genastropyskypointings(t_start,
                                t_stop,
                                dirs;
                                latitude_deg=0.0,
                                longitude_deg=0.0,
                                height_m=0.0)

    t_start_s = astropy_time[:Time](t_start, format="iso", scale="utc")
    t_stop_s = astropy_time[:Time](t_stop, format="iso", scale="utc")
    
    jd_range = range(t_start_s[:jd], stop=t_stop_s[:jd], length=size(dirs)[1])
    loc = astropy_coordinates[:EarthLocation](lon=longitude_deg,
                                              lat=latitude_deg,
                                              height=height_m)
    
    skydirs = Array{Float64}(undef, size(dirs))

    for (idx, time_jd) = enumerate(jd_range)
        Alt_rad = π/2 - dirs[idx, 1] 
        Az_rad = 2π - dirs[idx, 2]

        skycoord = astropy_coordinates[:SkyCoord](
            alt=Alt_rad*astropy_units[:rad],
            az=Az_rad*astropy_units[:rad],
            obstime=astropy_time[:Time](time_jd, format="jd"),
            frame="altaz",
            location=loc)
        
        skydirs[idx, 1] = deg2rad(skycoord[:transform_to]("icrs")[:ra][1]) 
        skydirs[idx, 2] = deg2rad(skycoord[:transform_to]("icrs")[:dec][1])

    end

    skydirs
end


"""
    genskypointings(t_start, t_stop, dirs; latitude_deg=0.0, longitude_deg=0.0, 
                    height_m=0.0)

Project a set of pointings in the sky.  

The parameter `t_start` and `t_start` must be two strings which tells the exact 
UTC date and time of the observation in "iso" format. The parameter `dirs` must 
be a N×2 array containing the observed directions expressed in colatitude and 
the longitude. The keywords `latitude_deg`, `longitude_deg` and `height_m` should
contain the latitude (in degrees, N is positive), the longitude (in degrees, 
counterclockwise is positive) and the height (in meters)  of the location where 
the observation is made.

Return a 2-tuple containing the observed directions projected in sky coordinates 
(a N×2 array containing the Right ascension and the Declination, in radians) at 
each time step. 
Directions are expressed in ICRS coordinates.

Example:
    `````julia
genskypointings("2019-01-01 00:00:00", "2020-01-25 14:52:10.05", dirs)
`````
"""
function genskypointings(t_start,
                         t_stop,
                         dirs;
                         latitude_deg=0.0,
                         longitude_deg=0.0,
                         height_m=0.0)

    t_start_s = astropy_time[:Time](t_start, format="iso", scale="utc")
    t_stop_s = astropy_time[:Time](t_stop, format="iso", scale="utc")
    
    jd_range = range(t_start_s[:jd], stop=t_stop_s[:jd], length=size(dirs)[1])
    loc = astropy_coordinates[:EarthLocation](lon=longitude_deg,
                                              lat=latitude_deg,
                                              height=height_m)
    
    skydirs = Array{Float64}(undef, size(dirs))

    for (idx, time_jd) = enumerate(jd_range)
        Alt_rad = π/2 - dirs[idx, 1] 
        Az_rad = 2π - dirs[idx, 2]

        LST = astropy_time[:Time](time_jd, format="jd",
                                  location=loc)[:sidereal_time]("mean")[1]
        Lat_rad = deg2rad(latitude_deg)
        Dec = asin(sin(Alt_rad) * sin(Lat_rad) + cos(Alt_rad) * cos(Lat_rad) *
                   cos(Az_rad))
        HourAngle = acos((sin(Alt_rad) - sin(Dec) * sin(Lat_rad)) /
                         (cos(Dec) * cos(Lat_rad)))
        
        if sin(Az_rad) < 0
            h = rad2deg(HourAngle) / 15
        else 
            h = (360 - rad2deg(HourAngle)) / 15
        end

        Ra = LST - h
        if Ra < 0 Ra += 24 end

        skydirs[idx, 1] = deg2rad(360 * Ra / 24)
        skydirs[idx, 2] = Dec
    end

    skydirs
end

