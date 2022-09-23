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

export configuration_angles_ST
export camtoground
export startrackerpointings, startrackerpointings!

"""
    configuration_angles(wheel1ang_0_rad,
                         wheel2ang_0_rad,
                         wheel3ang_0_rad,
                         forkang_rad,
                         omegaVAXang_rad,
                         zVAXang_rad,
                         rollang_rad,
                         panang_rad,
                         tiltang_rad)

Struct containing the configuration angles for the star tracker ie the angles describing
the non idealities in the telescope and in the camera reference frame(all of these parameters are ù
considered equal to 0 in an ideal case):

(rollang, panang, tiltang): Tait-Brian angles encoding the camera orientation in the telescope reference frame

(wheel1ang_0, wheel2ang_0, wheel3ang_0): these are the zero points angles for the three motors
                                         (respectively the boresight, the altitude and the ground
                                         motor) of the telescope

(forkang): describe the deviation of orthogonality between the H-AXIS and the V-AXIS

(omegaVAXang, zVAXang): wobble angles encoding the deviation of the V-AXIS from the local vertical;
                        zVAXang is the displacement from the V-AXIS ie the colatitude,
                        omegaVAXang is the azimuth of the ascending node defined as 
                        omegaVAXang = 90° + A * zVAXang
                        
All of these angles must be expressed in RADIANS and measured anticlockwise.
"""
Base.@kwdef struct configuration_angles_ST
    wheel1ang_0_rad :: Float64 = 0
    wheel2ang_0_rad :: Float64 = 0
    wheel3ang_0_rad :: Float64 = 0
    forkang_rad :: Float64 = 0
    omegaVAXang_rad :: Float64 = 0
    zVAXang_rad :: Float64 = 0
    rollang_rad :: Float64 = 0
    panang_rad :: Float64 = 0
    tiltang_rad :: Float64 = 0
end

"""
    camtoground(wheelanglesfn,
                time_s,
                config_st)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transformation for the star tracker, taking into account
all the non idealities: the position of the camera in the telescope 
reference frame and the configuration angles of the telescope.
The parameter `config_st` must be a configuration_angles_ST struct
containing the angles describing the non idealities of the telescope.
"""
function camtoground(wheelanglesfn, 
                     time_s,
                     config_st::configuration_angles_ST = configuration_angles_ST())
    (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

    qroll = qrotation_z(config_st.rollang_rad)
    qpan = qrotation_y(config_st.panang_rad)
    qtilt = qrotation_x(config_st.tiltang_rad)

    qwheel1 = qrotation_z(wheel1ang - config_st.wheel1ang_0_rad)
    qwheel2 = qrotation_y(wheel2ang - config_st.wheel2ang_0_rad)
    # The minus sign here takes into account the fact that the azimuth
    # motor requires positive angles to turn North into East
    qwheel3 = qrotation_z(-wheel3ang + config_st.wheel3ang_0_rad)

    qfork = qrotation_x(config_st.forkang_rad)
    qomegaVAX = qrotation_z(config_st.omegaVAXang_rad)
    qzVAX = qrotation_x(config_st.zVAXang_rad)    

    qomegaVAX * (qzVAX * (qwheel3 * (qfork * (qwheel2 * (qwheel1 * ( qtilt * (qpan * qroll)))))))
end

function startrackerpointings!(wheelanglesfn,
                               config_ang::configuration_angles_ST,
                               st_dir,
                               timerange_s,
                               dirs,
                               latitude_deg = TENERIFE_LATITUDE_DEG,
                               ground = false,
                               day_duration_s = 86400.0)

    if ground
        @assert size(dirs, 2) == 4
    else
        @assert size(dirs, 2) == 2
    end

    for (idx, time_s) = enumerate(timerange_s)
        # This converts the RDP into the MCS (ground reference frame)
        groundq = camtoground(wheelanglesfn, time_s, config_ang)
        # This converts the MCS into the celestial reference frame
        quat = groundtoearth(groundq, time_s, latitude_deg; day_duration_s = day_duration_s)
        θ, ϕ, _ = quat_to_angles(st_dir, Float64[1.0, 0.0, 0.0], quat)
        (dirs[idx, 1], dirs[idx, 2]) = (θ, ϕ)

        if ground
            # Re-run the transformation algorithm using the ground quaternion
            θ_ground, ϕ_ground, _ = quat_to_angles(st_dir, Float64[1.0, 0.0, 0.0], groundq)

            (dirs[idx, 3], dirs[idx, 4]) = (θ_ground, ϕ_ground)
        end
    end
end

function startrackerpointings(wheelanglesfn,
                              config_ang::configuration_angles_ST,
                              st_dir,
                              timerange_s;
                              latitude_deg = TENERIFE_LATITUDE_DEG,
                              ground = false,
                              day_duration_s = 86400.0)

    if ground
        dirs = Array{Float64}(undef, length(timerange_s), 4)
    else
        dirs = Array{Float64}(undef, length(timerange_s), 2)
    end



    dirs
end