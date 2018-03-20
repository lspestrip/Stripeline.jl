using Quaternions
import Healpix

export TENERIFE_LATITUDE_DEG, timetorotang, genpointings

TENERIFE_LATITUDE_DEG = 28.29


"""
    timetorotang(time, rpm)

Convert a time into a rotation angle, given the number of rotations per minute.
The time should be expressed in seconds. The return value is in radians.
`time` can either be a scalar or a vector.
"""

function timetorotang(time, rpm)
    if rpm == 0
        0.0
    else
        2 * π * time * (rpm / 60)
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

Example:
`````julia
genpointings([0, 0, 1], 0:0.1:1) do time_s
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    return (0.0, 20.0, timetorotang(time_s, 1))
end
`````
"""

function genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)
    
    dirs = Array{Float64}(length(timerange_s), 2)
    ψ = Array{Float64}(length(timerange_s))

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
        poldir = rotmatr * [1; 0; 0]
        
        # The North for a vector v is just dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(vector[1], vector[2], vector[3])
        dirs[idx, :] = [θ, ϕ]
        northdir = [cos(θ) * cos(ϕ); cos(θ) * sin(ϕ); sin(θ)]
        
        cosψ = clamp(dot(northdir, poldir), -1, 1)
        crosspr = northdir × poldir
        sinψ = clamp(dot(crosspr, crosspr), -1, 1)
        ψ[idx] = atan2(cosψ, sinψ) * sign(dot(crosspr, vector))
    end
    
    (dirs, ψ)
end