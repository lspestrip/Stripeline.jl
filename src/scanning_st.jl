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

include("quaternions.jl")


"""
    configuration_angles(wheel1ang_0,
                         wheel2ang_0,
                         wheel3ang_0,
                         forkang,
                         omegaVAXang,
                         zVAXang,
                         rollang,
                         panang,
                         tiltang)

Struct containing the configuration angles for the star tracker ie the angles describing
the non idealities in the telescope and in the camera reference frame(all of these parameters are ù
considered equal to 0 in an ideal case):

(rollang, panang, tiltang): 

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
Base.@kwdef struct configuration_angles_ST
    wheel1ang_0 :: Float64 = 0
    wheel2ang_0 :: Float64 = 0
    wheel3ang_0 :: Float64 = 0
    forkang :: Float64 = 0
    omegaVAXang :: Float64 = 0
    zVAXang :: Float64 = 0
    rollang :: Float64 = 0
    panang :: Float64 = 0
    tiltang :: Float64 = 0
end