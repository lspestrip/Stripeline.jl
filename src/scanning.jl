# -*- encoding: utf-8 -*-
#
# MIT License
#
# Copyright (c) 2016 Maurizio Tomasi
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
# # TERMINOLOGY
#
# The code in this file tries to follow the conventions laid in
# document LSPE-STRIP-SP-017 ("STRIP Telescope Coordinate Systems
# Specification") as much as possible. As the reader might not be
# acquainted with it, here is a short summary of the relevant
# terminology:
#
# +-----+-----------------------------+------------------------------------------+
# | SCS | Site Coordinate System      | Centered on the GPS station              |
# +-----+-----------------------------+------------------------------------------+
# | MCS | Mount Coordinate System     | Centered at the basis of the telescope   |
# +-------+---------------------------+------------------------------------------+
# | ACS | Azimuth Coordinate System   | Moveable, centered in the middle of      |
# |     |                             | ground motor. Positive angles go from    |
# |     |                             | North to *East* (beware!)                |
# +-----+-----------------------------+------------------------------------------+
# | ECS | Elevation Coordinate System | Moveable, the Y axis is aligned with     |
# |     |                             | the motor axis                           |
# +-----+-----------------------------+------------------------------------------+
# | BCS | Boresight Coordinate System | It should have been moveable, but the    |
# |     |                             | current design does not allow it (alas!) |
# +-----+-----------------------------+------------------------------------------+
# | RDP | Reference Detector Plane    | Centered in the focus of the telescope.  |
# |     |                             | All the beam vectors are specified in    |
# |     |                             | this frame.                              |
# +-----+-----------------------------+------------------------------------------+
#
# We always use a right-handed coordinate system:
#
#                 ^ Z axis (aligned with "up" direction, e.g., Zenith)
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#          Origin O----------------> Y axis
#                /
#               /
#              /
#             /
#            /
#           v  X axis (aligned with North, when applicable)
#

import Healpix
import StaticArrays
import LinearAlgebra: ×, dot
import AstroLib
import Dates
include("quaternions.jl")

export TENERIFE_LATITUDE_DEG, TENERIFE_LONGITUDE_DEG, TENERIFE_HEIGHT_M
export TelescopeAngles, CameraAngles
export directiontoangles
export timetorotang, camtotelescope, telescopetoground, groundtoearth
export genpointings!, genpointings, northdir, eastdir, polarizationangle

"Latitude of the LSPE/Strip site in Tenerife, in degrees"
const TENERIFE_LATITUDE_DEG = 28.30026

"Longitude of the LSPE/Strip site in Tenerife, in degrees"
const TENERIFE_LONGITUDE_DEG = -16.51012

"Height of the LSPE/Strip site in Tenerife, in meters"
const TENERIFE_HEIGHT_M = 2390

"""
    TelescopeAngles(
        wheel1ang_0_rad :: Float64 = 0.0,
        wheel2ang_0_rad :: Float64 = 0.0,
        wheel3ang_0_rad :: Float64 = 0.0,
        forkang_rad :: Float64 = 0.0,
        omegaVAXang_rad :: Float64 = 0.0,
        zVAXang_rad :: Float64 = 0.0  
    )

Struct containing the configuration angles for the telescope i.e. the angles describing
the non idealities in the telescope (all of these parameters are considered equal to 0 in 
an ideal telescope):

(`wheel1ang_0_rad`, `wheel2ang_0_rad`, `wheel3ang_0_rad`): these are the zero points angles for the three motors
                                                           (respectively the boresight, the altitude and the ground
                                                           motor)

(`forkang_rad`): describe the deviation of orthogonality between the H-AXIS and the V-AXIS

(`omegaVAXang_rad`, `zVAXang_rad`): wobble angles encoding the deviation of the V-AXIS from the local vertical;
                                    zVAXang is the displacement from the V-AXIS,
                                    omegaVAXang is the azimuth of the ascending node.

See the documentation for a graphical rapresentation of each angles.

All of these angles must be expressed in RADIANS and measured anticlockwise.
"""
Base.@kwdef struct TelescopeAngles
    wheel1ang_0_rad :: Float64 = 0.0
    wheel2ang_0_rad :: Float64 = 0.0
    wheel3ang_0_rad :: Float64 = 0.0
    forkang_rad :: Float64 = 0.0
    omegaVAXang_rad :: Float64 = 0.0
    zVAXang_rad :: Float64 = 0.0
end

"""
    CameraAngles(
        panang_rad :: Float64 = 0.0,
        tiltang_rad :: Float64 = 0.0,
        rollang_rad :: Float64 = 0.0
    )

Struct encoding the Tait-Bryan angles of a single camera/detector.
Camera angles represent the orientation of the detector relative 
to the telescope.

- `panang_rad` encode a rotation around x-axis
- `tiltang_rad` encode a rotation around y-axis
- `rollang_rad` encode a rotation around z-axis

See [`camtotelescope`](@ref) to understand how these angles are used to
transform the pointing direction.

All of these angles must be expressed in RADIANS and measured anticlockwise.
"""
Base.@kwdef struct CameraAngles
    panang_rad :: Float64 = 0.0
    tiltang_rad :: Float64 = 0.0
    rollang_rad :: Float64 = 0.0
end

"""
    directiontoangles(dir)

This function convert a pointing direction vector into Tait-Brian angles
used in [`CameraAngles`](@ref).

The convention used to rotate camera is R_x*R_y*R_z (see [`camtotelescope`](@ref)).

The input vector `dir` must be normalized.

This function is used internally by [`genpointings`](@ref) to
maintain backwards compatibility.
"""
function directiontoangles(dir)
    #y-axis rotation angle
    tiltang = asin(dir[1])
    #x-axis rotation angle
    panang = - asin(dir[2]/cos(tiltang))
    #z-azis rotation angle
    rollang = 0.0
    (panang, tiltang, rollang)
end

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
    camtotelescope(cam_ang::Camera_angles)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transform of the detector direction into the focal plane 
reference frame.

This function is used internally in [`genpointings`](@ref) as
a part of the rotations chain.
"""
function camtotelescope(cam_ang::CameraAngles)
    qroll = qrotation_z(cam_ang.rollang_rad)
    qtilt = qrotation_y(cam_ang.tiltang_rad)
    qpan = qrotation_x(cam_ang.panang_rad)
    qpan * (qtilt * qroll)
end


function telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang::Nothing = nothing)

    qwheel1 = qrotation_z(wheel1ang)
    qwheel2 = qrotation_y(wheel2ang)

    # The minus sign here takes into account the fact that the azimuth
    # motor requires positive angles to turn North into East
    qwheel3 = qrotation_z(-wheel3ang)

    qwheel3 * (qwheel2 * qwheel1)
end

function telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang::TelescopeAngles)

    qwheel1 = qrotation_z(wheel1ang - telescope_ang.wheel1ang_0_rad)
    qwheel2 = qrotation_y(wheel2ang - telescope_ang.wheel2ang_0_rad)
    qwheel3 = qrotation_z(-wheel3ang + telescope_ang.wheel3ang_0_rad)

    qfork = qrotation_x(telescope_ang.forkang_rad)
    qomegaVAX = qrotation_z(telescope_ang.omegaVAXang_rad)
    qzVAX = qrotation_x(telescope_ang.zVAXang_rad)    

    qomegaVAX * (qzVAX * (qwheel3 * (qfork * (qwheel2 * qwheel1))))
end 

@doc raw"""
    telescopetoground(wheelanglesfn, time_s, telescope_ang::Nothing = nothing)
    telescopetoground(wheelanglesfn, time_s, telescope_ang::TelescopeAngles)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transform from the focal plane to the ground of the
telescope. The parameter `wheelanglesfn` must be a function which
takes as input a time, `time_s`, in seconds, and it must return a
3-tuple containing the angles of the following motors, measured in
**radians**:

1. The boresight motor (rotation around the ``z`` axis,
   counterclockwise)

2. The altitude motor (rotation around the ``y`` axis,
   counterclockwise)

3. The ground motor (rotation around the ``z`` axis, **clockwise**:
   N→E→S→W)

The parameter telescope_ang must be a `TelescopeAngles` struct containing the angles 
describing the non idealities of the telescope. If `nothing` is passed the function
calculate the quaternion associated with the ideal case i.e. like all the TelescopeAngles
are zero.

# Example

`````julia
telescopetoground(3600.0) do
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````
"""
telescopetoground

"""
    groundtoearth(groundq, time_s, latitude_deg; day_duration_s=86400.0)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transformation from the ground of the telescope to the
Equatorial coordinate system. The parameter `groundq` must be a
quaternion describing the coordinate transformation from the focal
plane of the telescope to the ground. The parameter `time_s` must be a
time in seconds, and `latitude_deg` is the latitude (in degrees, N is
positive) of the location where the observation is made.

The keyword `day_duration_s` specifies the length of a sidereal day in
seconds.
"""
function groundtoearth(groundq, time_s, latitude_deg; day_duration_s = 86400.0)
    locq = qrotation_y(deg2rad(90 - latitude_deg))
    earthq = qrotation_z(2π * time_s / day_duration_s)

    earthq * (locq * groundq)
end


"""
    vector2equatorial(dir, jd, latitude_deg, longitude_deg, height_m;
                      prec=true, nut=true, aber=true)

Transform the Healpix coordinates of a vector into equatorial coordinates.
The parameter `vector` is 3D vector, `jd` is the julian date. The
paramters `latitude_deg`, `longitude_deg` and `height_m` should contain the
latitude (in degrees, N is positive), the longitude (in degrees, counterclockwise
is positive) and the height (in meters) of the location where the observation is
made. The parameters `prec`, `nut` and `aber` allow to consider secondary effects
respectively of precession, nutation and aberration.
"""
function vector2equatorial(vector, jd, latitude_deg, longitude_deg, height_m,
                           prec, nut, aber, ref)
    (θ, ϕ) = Healpix.vec2ang(vector...)
    alt_rad = π / 2 - θ
    az_rad = 2π - ϕ

    ra_deg, dec_deg, _ = AstroLib.hor2eq(rad2deg(alt_rad),
                                         rad2deg(az_rad),
                                         jd,
                                         latitude_deg,
                                         longitude_deg,
                                         height_m,
                                         precession = prec,
                                         nutate = nut,
                                         aberration = aber,
                                         refract = ref)
    (deg2rad(dec_deg), deg2rad(ra_deg))
end


"""
    polarizationangle(northdir, eastdir, poldir)

Calculate the polarization angle projected in the sky in IAU
conventions. The parameters `northdir` and `eastdir` must be versors
that point towards the North and East, respectively; `poldir` must be
a versor that identify the polarization direction projected in the
sky. The return value is in radians, and it is zero if the
polarization angles points toward East, π/2 if it points toward North,
etc.

# Examples
```jldoctest
julia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 1, 0])
0.0
julia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 0, 1]) |> rad2deg
90.0
julia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 0, -1]) |> rad2deg
-90.0
```
"""
function polarizationangle(northdir, eastdir, poldir)
    cosψ = clamp(dot(eastdir, poldir), -1, 1)
    sinψ = clamp(dot(northdir, poldir), -1, 1)
    atan(sinψ, cosψ)
end


northdir(θ, ϕ) = StaticArrays.@SArray [-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ)]
eastdir(θ, ϕ) = StaticArrays.@SArray [-sin(ϕ), cos(ϕ), 0.0]

"""
    northdir(θ, ϕ)
    eastdir(θ, ϕ)

Compute the North/East versor for a vector. The North for a vector v
is along the direction -dv/dθ, as θ is the colatitude and moves along
the meridian, and the East is along dv/dϕ.

# Examples
```jldoctest
julia> northdir(π/2, 0) ≈ [0, 0, 1]
true
julia> eastdir(π/2, 0) ≈ [0, 1, 0]
true
```
"""
northdir, eastdir


"""
    quat_to_angles(boreaxis, polaxis, quat)

Transform the boresight direction `dir` and the polarization direction
`poldir` according to quaternion `quat`. Return the 3-tuple (θ, ϕ, ψ)
representing the colatitude, longitude, and polarization angle
(calculated northward).

*Warning:* the definition of polarization angle does not work if the
observer is at the North or South Pole of the coordinate system. This
means that the polarization angle is not valid if the quaternion is
expressed in the ground reference system.

This function is used internally by [`genpointings`](@ref).
"""
function quat_to_angles(boreaxis::Array, polaxis::Array, quat::Quaternion)
    rotmatr = rotationmatrix_normalized(quat)
    boresight = rotmatr * boreaxis
    poldir = rotmatr * polaxis

    (θ, ϕ) = Healpix.vec2ang(boresight[1], boresight[2], boresight[3])

    north = northdir(θ, ϕ)
    east = eastdir(θ, ϕ)
    (θ, ϕ, polarizationangle(north, east, poldir))
end

"""
    quat_to_angles(quat)

Transform the boresight direction `[0.0,0.0,1.0]` and the polarization direction
`[1.0,0.0,0.0]` according to quaternion `quat`. Return the 3-tuple (θ, ϕ, ψ)
representing the colatitude, longitude, and polarization angle
(calculated northward).

*Warning:* the definition of polarization angle does not work if the
observer is at the North or South Pole of the coordinate system. This
means that the polarization angle is not valid if the quaternion is
expressed in the ground reference system.

This function is used internally by [`genpointings`](@ref).
"""
function quat_to_angles(quat::Quaternion)
    boresight = rotate_zaxis(quat)
    poldir = rotate_xaxis(quat)

    (θ, ϕ) = Healpix.vec2ang(boresight[1], boresight[2], boresight[3])

    north = northdir(θ, ϕ)
    east = eastdir(θ, ϕ)
    (θ, ϕ, polarizationangle(north, east, poldir))
end

# Old version of genpointings that accept beam_dir as a versor (Array) representing the pointing direction
function genpointings!(wheelanglesfn,
                       beam_dir,
                       timerange_s,
                       dirs,
                       psi;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       ground = false,
                       day_duration_s = 86400.0,
                       telescope_ang::Union{TelescopeAngles, Nothing} = nothing)

    if ground
        @assert size(dirs, 2) == 4
        @assert size(psi, 2) == 2
    else
        @assert size(dirs, 2) == 2
        @assert size(psi, 2) == 1
    end
    @assert size(dirs, 1) == size(psi, 1)

    for (idx, time_s) = enumerate(timerange_s)
        
        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

        # This converts the RDP into the MCS (ground reference frame)
        groundq = telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang)
        # This converts the MCS into the celestial reference frame
        quat = groundtoearth(groundq, time_s, latitude_deg; day_duration_s = day_duration_s)
            
        θ, ϕ, curpsi = quat_to_angles(beam_dir, polaxis, quat)
        (dirs[idx, 1], dirs[idx, 2]) = (θ, ϕ)

        if ground
            # Re-run the transformation algorithm using the ground quaternion
            θ_ground, ϕ_ground, psi_ground = quat_to_angles(beam_dir, polaxis, groundq)

            (dirs[idx, 3], dirs[idx, 4]) = (θ_ground, ϕ_ground)
            (psi[idx, 1], psi[idx, 2]) = (curpsi, psi_ground)
        else
            psi[idx] = curpsi
        end
    end
end

# New version of genpointings using CameraAngles
function genpointings!(wheelanglesfn,
                       beam_dir::CameraAngles,
                       timerange_s,
                       dirs,
                       psi;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       ground = false,
                       day_duration_s = 86400.0,
                       telescope_ang::Union{TelescopeAngles, Nothing} = nothing)

    if ground
        @assert size(dirs, 2) == 4
        @assert size(psi, 2) == 2
    else
        @assert size(dirs, 2) == 2
        @assert size(psi, 2) == 1
    end
    @assert size(dirs, 1) == size(psi, 1)

    camtotel_quat = camtotelescope(beam_dir)

    for (idx, time_s) = enumerate(timerange_s)

        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

        # This converts the RDP into the MCS (ground reference frame)
        groundq = telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang) * camtotel_quat
        # This converts the MCS into the celestial reference frame
        quat = groundtoearth(groundq, time_s, latitude_deg; day_duration_s = day_duration_s)
            
        θ, ϕ, curpsi = quat_to_angles(quat)
        (dirs[idx, 1], dirs[idx, 2]) = (θ, ϕ)

        if ground
            # Re-run the transformation algorithm using the ground quaternion
            θ_ground, ϕ_ground, psi_ground = quat_to_angles(groundq)

            (dirs[idx, 3], dirs[idx, 4]) = (θ_ground, ϕ_ground)
            (psi[idx, 1], psi[idx, 2]) = (curpsi, psi_ground)
        else
            psi[idx] = curpsi
        end
    end
end

function genpointings(wheelanglesfn,
                       beam_dir,
                       timerange_s;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       ground = false,
                       day_duration_s = 86400.0,
                       telescope_ang::Union{TelescopeAngles, Nothing} = nothing)

    if ground
        dirs = Array{Float64}(undef, length(timerange_s), 4)
        psi = Array{Float64}(undef, length(timerange_s), 2)
    else
        dirs = Array{Float64}(undef, length(timerange_s), 2)
        psi = Array{Float64}(undef, length(timerange_s))
    end

    genpointings!(
        wheelanglesfn,
        beam_dir,
        timerange_s,
        dirs,
        psi;
        polaxis = polaxis,
        latitude_deg = latitude_deg,
        ground = ground,
        day_duration_s = day_duration_s,
        telescope_ang = telescope_ang
    )

    (dirs, psi)
end

# Older version of genpointings using beam_dir::Array representing the pointing direction of the detector
function genpointings!(wheelanglesfn,
                       beam_dir,
                       timerange_s,
                       t_start::Dates.DateTime,
                       skydirs,
                       skypsi;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       longitude_deg = TENERIFE_LONGITUDE_DEG,
                       height_m = TENERIFE_HEIGHT_M,
                       precession = true,
                       nutation = true,
                       aberration = true,
                       refraction = true,
                       telescope_ang::Union{TelescopeAngles, Nothing} = nothing)

    @assert size(skydirs, 1) == size(skypsi, 1)
    @assert size(skydirs, 2) == 2
    @assert size(skypsi, 2) == 1

    for (idx, time_s) = enumerate(timerange_s)
        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)
        groundq = telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang)
        rotmatr = rotationmatrix_normalized(groundq)
        vector = rotmatr * beam_dir

        jd = AstroLib.jdcnv(t_start + Dates.Nanosecond(round(Int64, time_s * 1e9)))
        Dec_rad, Ra_rad = vector2equatorial(vector,
                                            jd,
                                            latitude_deg,
                                            longitude_deg,
                                            height_m,
                                            precession,
                                            nutation,
                                            aberration,
                                            refraction)

        poldir = rotmatr * polaxis
        north = northdir(π / 2 - Dec_rad, Ra_rad)
        east = eastdir(π / 2 - Dec_rad, Ra_rad)

        skydirs[idx, 1] = π/2 - Dec_rad
        skydirs[idx, 2] = Ra_rad
        skypsi[idx] = polarizationangle(north, east, poldir)
    end
end

# New version using CameraAngles
function genpointings!(wheelanglesfn,
                       beam_dir::CameraAngles,
                       timerange_s,
                       t_start::Dates.DateTime,
                       skydirs,
                       skypsi;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       longitude_deg = TENERIFE_LONGITUDE_DEG,
                       height_m = TENERIFE_HEIGHT_M,
                       precession = true,
                       nutation = true,
                       aberration = true,
                       refraction = true,
                       telescope_ang::Union{TelescopeAngles, Nothing} = nothing)

    @assert size(skydirs, 1) == size(skypsi, 1)
    @assert size(skydirs, 2) == 2
    @assert size(skypsi, 2) == 1

    camtotel_quat = camtotelescope(beam_dir)

    for (idx, time_s) = enumerate(timerange_s)
        (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)
        groundq = telescopetoground(wheel1ang, wheel2ang, wheel3ang, telescope_ang) * camtotel_quat
        vector = rotate_zaxis(groundq)

        jd = AstroLib.jdcnv(t_start + Dates.Nanosecond(round(Int64, time_s * 1e9)))
        Dec_rad, Ra_rad = vector2equatorial(vector,
                                            jd,
                                            latitude_deg,
                                            longitude_deg,
                                            height_m,
                                            precession,
                                            nutation,
                                            aberration,
                                            refraction)

        poldir = rotate_xaxis(groundq)
        north = northdir(π / 2 - Dec_rad, Ra_rad)
        east = eastdir(π / 2 - Dec_rad, Ra_rad)

        skydirs[idx, 1] = π/2 - Dec_rad
        skydirs[idx, 2] = Ra_rad
        skypsi[idx] = polarizationangle(north, east, poldir)
    end
end

function genpointings!(wheelanglesfn,
                       beam_dir,
                       timerange_s,
                       t_start::Number,
                       skydirs,
                       skypsi;
                       kwargs...)
    genpointings!(
        wheelanglesfn,
        beam_dir,
        timerange_s,
        AstroLib.daycnv(t_start)::Dates.DateTime,
        skydirs,
        skypsi,
        kwargs...)
end

function genpointings(wheelanglesfn,
                      beam_dir,
                      timerange_s,
                      t_start::Dates.DateTime;
                      kwargs...)

    skydirs = Array{Float64}(undef, length(timerange_s), 2)
    skypsi = Array{Float64}(undef, length(timerange_s))

    genpointings!(
        wheelanglesfn,
        beam_dir,
        timerange_s,
        t_start,
        skydirs,
        skypsi;
        kwargs...)

    (skydirs, skypsi)
end

function genpointings(wheelanglesfn,
                       beam_dir,
                       timerange_s,
                       t_start::Number;
                       kwargs...)
    genpointings(
        wheelanglesfn,
        beam_dir,
        timerange_s,
        AstroLib.daycnv(t_start)::Dates.DateTime;
        kwargs...)
end


@doc raw"""
    genpointings!(wheelanglesfn, beam_dir::CameraAngles, 
                  timerange_s, 
                  dirs, psi;
                  polaxis = Float64[1.0, 0.0, 0.0],
                  latitude_deg = TENERIFE_LATITUDE_DEG,
                  ground = false,
                  config_ang::Union{ConfigAngles, Nothing} = nothing)
    genpointings(wheelanglesfn, beam_dir::CameraAngles, 
                 timerange_s;
                 polaxis = Float64[1.0, 0.0, 0.0],
                 latitude_deg = TENERIFE_LATITUDE_DEG,
                 ground = false,
                 config_ang::Union{ConfigAngles, Nothing} = nothing)
    genpointings!(wheelanglesfn, beam_dir::CameraAngles, 
                  timerange_s, t_start, dirs, psi;
                  polaxis = Float64[1.0, 0.0, 0.0],
                  latitude_deg = TENERIFE_LATITUDE_DEG,
                  longitude_deg = TENERIFE_LONGITUDE_DEG,
                  height_m = TENERIFE_HEIGHT_M,
                  precession = true,
                  nutation = true,
                  aberration = true,
                  refraction = true,
                  config_ang::Union{ConfigAngles, Nothing} = nothing)
    genpointings(wheelanglesfn, beam_dir::CameraAngles, 
                 timerange_s, t_start;
                 polaxis=Float64[1.0, 0.0, 0.0],
                 latitude_deg=TENERIFE_LATITUDE_DEG,
                 longitude_deg=TENERIFE_LONGITUDE_DEG,
                 height_m=TENERIFE_HEIGHT_M,
                 precession = true,
                 nutation = true,
                 aberration = true,
                 refraction = true,
                 config_ang::Union{ConfigAngles, Nothing} = nothing)

Generate a set of pointing directions for a STRIP detector. Each
function is provided in two flavours: the ones ending with `!` save
the results in the last two parameters `dirs` and `psi`, while the
others automatically allocate memory and return their results as a
pair `(dirs, psi)`.

The parameter `wheelanglesfn` must be a function which takes as input
a time in seconds and returns a 3-tuple containing the angles (in
radians) of the three motors:

1. The boresight motor (rotation around the ``z`` axis,
   counterclockwise)

2. The altitude motor (rotation around the ``y`` axis,
   counterclockwise)

3. The ground motor (rotation around the ``z`` axis, **clockwise**:
   N→E→S→W)

The meaning of the parameters/keywords is the following:

- `beam_dir` is [`CameraAngles`](@ref) object that specifies the three
  Tait-Bryan angles for the detector.
  
  N.B. in a previous version this parameter specifies the normalized 
  pointing direction of the mean (the boresight is [0, 0, 1]). This is
  still possible, due to an internal conversion if we use an `Array` as beam_dir, 
  but will no longer be supported as default behavior.

- `timerange_s` is an enumerable type that specifies at which times
  (in seconds) pointings must be computed.

- `t_start` is a DateTime object or a Julian date, which specifies the
  UTC date and time when the observation starts

- `latitude_deg` is the latitude of the location where the observation
  is made (in degrees, North is positive). The default value is
  [`TENERIFE_LATITUDE_DEG`](@ref).

- `longitude_deg` is the longitude of the location where the
  observation is made. The default value is
  [`TENERIFE_LONGITUDE_DEG`](@ref).

- `height_m` is the elevation of the location where the observation is
  made (in meters). The default value is [`TENERIFE_HEIGHT_M`](@ref).

- If `ground` is `true`, the function will return a 4-tuple containing
  the colatitude and longitude measured in Equatorial coordinates
  (columns 1 and 2) and in ground coordinates (columns 3 and 4). If
  `ground` is `false`, only the Equatorial coordinates are computed.

- `polaxis` is the polarization axis; it must be normalized.

- `precession`: if `true` (the default), the Earth's precession is
  taken into account.

- `nutation`: if `true` (the default), the Earth's nutation is taken
  into account.

- `aberration`: if `true` (the default), stellar aberration is taken
  into account.

- `refraction`: if `true`, refraction corrections are taken into account.
  As these corrections are only valid for optical wavelengths, the
  default is `false`.

- `telescope_ang`: specifies the configuration angles for the telescope 
  (see [`TelescopeAngles`](@ref) for more details). This is used
  internally by [`telescopetoground`](@ref); if nothing is passes then the version 
  of telescopetoground for an ideal telescope will be used.

# Return values

If `t_start` is not provided, the function `genpointings` returns a
2-tuple containing the sky directions (a N×2 array containing
declination and right ascension, in Equatorial coordinates) and the
polarization angle for each time step. The function `genpointings!` sets
the values in the last two parameters `dirs` and `psi`.

If `t_start` is provided, the function `genpointings` returns a
2-tuple (4-tuple) containing the directions (a N×2 or Nx4 array
containing the declination and the right ascension) and the polarization
angles at each time step; `genpointings!` works as above.


# Examples

Here is an example using the form without `t_start`:

`````julia
dir, psi = genpointings(CameraAngles(), 0:0.1:1) do time_s
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````

And here is an example using `t_start`; unlike the previous example,
we use a lambda function instead of the `do...end` notation.

`````julia
import Dates

dirs, psi = genpointings(time_s -> (0, deg2rad(20),
                                    timetorotang(time_s, 1)),
                         CameraAngles(),
                         0:0.1:1,
                         Dates.DateTime(2022, 01, 01, 0, 0, 0),
                         latitude_deg=10.0,
                         longitude_deg=20.0,
                         height_m = 1000) do time_s
`````
"""
genpointings
