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

export TENERIFE_LATITUDE_DEG, TENERIFE_LONGITUDE_DEG, TENERIFE_HEIGHT_M
export configuration_angles, configuration_angles_ST
export timetorotang, telescopetoground, groundtoearth
export genpointings!, genpointings, northdir, eastdir, polarizationangle

"Latitude of the LSPE/Strip site in Tenerife, in degrees"
const TENERIFE_LATITUDE_DEG = 28.30026

"Longitude of the LSPE/Strip site in Tenerife, in degrees"
const TENERIFE_LONGITUDE_DEG = -16.51012

"Height of the LSPE/Strip site in Tenerife, in meters"
const TENERIFE_HEIGHT_M = 2390

include("quaternions.jl")

"""
    configuration_angles(wheel1ang_0,
                         wheel2ang_0,
                         wheel3ang_0,
                         forkang,
                         omegaVAXang,
                         zVAXang)

Struct containing the configuration angles for the telescope ie the angles describing
the non idealities in the telescope (all of these parameters are considered equal to 0 in 
an ideal telescope):

(wheel1ang_0, wheel2ang_0, wheel3ang_0): these are the zero points angles for the three motors
                                         (respectively the boresight, the altitude and the ground
                                         motor)

(forkang): describe the deviation of orthogonality between the H-AXIS and the V-AXIS

(omegaVAXang, zVAXang): wobble angles encoding the deviation of the V-AXIS from the local vertical;
                        zVAXang is the displacement from the V-AXIS ie the colatitude,
                        omegaVAXang is the azimuth of the ascending node defined as 
                        omegaVAXang = 90° + A * zVAXang
                        
All of these angles must be expressed in RADIANS.
"""
Base.@kwdef struct configuration_angles
    wheel1ang_0 :: Float64 = 0
    wheel2ang_0 :: Float64 = 0
    wheel3ang_0 :: Float64 = 0
    forkang :: Float64 = 0
    omegaVAXang :: Float64 = 0
    zVAXang :: Float64 = 0
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


@doc raw"""
    telescopetoground(wheelanglesfn, time_s)

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
function telescopetoground(wheelanglesfn, time_s)
    (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

    qwheel1 = qrotation_z(wheel1ang)
    qwheel2 = qrotation_y(wheel2ang)

    # The minus sign here takes into account the fact that the azimuth
    # motor requires positive angles to turn North into East
    qwheel3 = qrotation_z(-wheel3ang)

    qwheel3 * (qwheel2 * qwheel1)
end


"""
    telescopetoground(wheelanglesfn, time_s, config_ang)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transform from the focal plane to the ground of the
telescope. 
The parameter `config_ang` must be a configuration_angles struct
containing the angles describing the non idealities of the telescope.

"""
function telescopetoground(wheelanglesfn, time_s, config_ang::configuration_angles = configuration_angles())
    (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

    qwheel1 = qrotation_z(wheel1ang - config_ang.wheel1ang_0)
    qwheel2 = qrotation_y(wheel2ang - config_ang.wheel2ang_0)

    # The minus sign here takes into account the fact that the azimuth
    # motor requires positive angles to turn North into East
    qwheel3 = qrotation_z(-wheel3ang + config_ang.wheel3ang_0)

    qfork = qrotation_x(config_ang.forkang)
    qomegaVAX = qrotation_z(config_ang.omegaVAXang)
    qzVAX = qrotation_x(config_ang.zVAXang)    

    qomegaVAX * (qzVAX * (qwheel3 * (qfork * (qwheel2 * qwheel1))))
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
function quat_to_angles(boreaxis, polaxis, quat)
    rotmatr = rotationmatrix_normalized(quat)
    boresight = rotmatr * boreaxis
    poldir = rotmatr * polaxis

    (θ, ϕ) = Healpix.vec2ang(boresight...)

    north = northdir(θ, ϕ)
    east = eastdir(θ, ϕ)
    (θ, ϕ, polarizationangle(north, east, poldir))
end


function genpointings!(wheelanglesfn,
                       beam_dir,
                       timerange_s,
                       dirs,
                       psi;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       ground = false,
                       day_duration_s = 86400.0)

    if ground
        @assert size(dirs, 2) == 4
        @assert size(psi, 2) == 2
    else
        @assert size(dirs, 2) == 2
        @assert size(psi, 2) == 1
    end
    @assert size(dirs, 1) == size(psi, 1)

    for (idx, time_s) = enumerate(timerange_s)

        # This converts the RDP into the MCS (ground reference frame)
        groundq = telescopetoground(wheelanglesfn, time_s)
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

function genpointings(wheelanglesfn,
                       beam_dir,
                       timerange_s;
                       polaxis = Float64[1.0, 0.0, 0.0],
                       latitude_deg = TENERIFE_LATITUDE_DEG,
                       ground = false,
                       day_duration_s = 86400.0)

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
    )

    (dirs, psi)
end

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
                       refraction = true)

    @assert size(skydirs, 1) == size(skypsi, 1)
    @assert size(skydirs, 2) == 2
    @assert size(skypsi, 2) == 1

    for (idx, time_s) = enumerate(timerange_s)
        groundq = telescopetoground(wheelanglesfn, time_s)
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

        skydirs[idx, 1] = Dec_rad
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
    genpointings!(wheelanglesfn, beam_dir, timerange_s, dirs, psi;
                  polaxis = Float64[1.0, 0.0, 0.0],
                  latitude_deg = TENERIFE_LATITUDE_DEG,
                  ground = false)
    genpointings(wheelanglesfn, beam_dir, timerange_s;
                 polaxis = Float64[1.0, 0.0, 0.0],
                 latitude_deg = TENERIFE_LATITUDE_DEG,
                 ground = false)
    genpointings!(wheelanglesfn, beam_dir, timerange_s, t_start, dirs, psi;
                  polaxis = Float64[1.0, 0.0, 0.0],
                  latitude_deg = TENERIFE_LATITUDE_DEG,
                  longitude_deg = TENERIFE_LONGITUDE_DEG,
                  height_m = TENERIFE_HEIGHT_M,
                  precession = true,
                  nutation = true,
                  aberration = true,
                  refraction = true)
    genpointings(wheelanglesfn, beam_dir, timerange_s, t_start;
                 polaxis=Float64[1.0, 0.0, 0.0],
                 latitude_deg=TENERIFE_LATITUDE_DEG,
                 longitude_deg=TENERIFE_LONGITUDE_DEG,
                 height_m=TENERIFE_HEIGHT_M,
                 precession = true,
                 nutation = true,
                 aberration = true,
                 refraction = true)

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

- `beam_dir` specifies the pointing direction of the mean (the boresight is
  [0, 0, 1]). It must be normalized.

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
dir, psi = genpointings([0, 0, 1], 0:0.1:1) do time_s
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
                         [0, 0, 1],
                         0:0.1:1,
                         Dates.DateTime(2022, 01, 01, 0, 0, 0),
                         latitude_deg=10.0,
                         longitude_deg=20.0,
                         height_m = 1000) do time_s
`````
"""
genpointings
