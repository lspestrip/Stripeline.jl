```@meta
DocTestSetup = quote
    using Stripeline
end
```

# [Pointing Reconstruction Model (PRM)](@id prm)

Starting from the scanning strategy described in [Scanning strategy](@ref scanning_strategy) we
can improve the model taking into account the non idealities of the system, parametrized by the so-called _configuration angles_.

The aim of the PRM is to calculate the pointing direction of a generic camera mounted on the STRIP
telescope, given the control angles (encoding the positions of the azimuth and altitudes motors)
as a function of time and the configuration angles describing the geometry of the telescope.

## Geometry of the telescope

First, to understand all the control angles lets analyze a model of the telescope and its 
alt-az mount. A simplified model is reported in the following figure.

```@raw html
<figure>
    <img src="../assets/prm/telescope_model.png" width="60%"/>
    <figcaption>Fig.1 Telescope model</figcaption>
</figure>
```

A basement holds a vertical axis (V-AXIS) allowing the Azimuth rotation, a fork mounted on the
top of the V-AXIS holds the horizontal axis (H-AXIS) which allows the Altitude rotation. In an
ideal case V-AXIS and H-AXIS are perpendicular and V-AXIS is aligned with the local
(topocentric) zenith.

The rotation around the V-AXIS is performed by the _ground_ motor, while the rotation around
H-AXIS by the _altitude_ motor. There is also a _boresight_ motor that will be kept constant
and equal to zero. The position of each motor is encoded by three angles (_wheel1ang_, _wheel2ang_, _wheel3ang_)
that describe, respectively, the boresight, the altitude and the ground motor rotation.

## Configuration angles

To describe the non idealities of the telescope we need six angles, each of which represent a rotation around
specific ($\hat{e}_x$, $\hat{e}_y$, $\hat{e}_z$) coordinate axis:

1. `wheel1ang_0`, `wheel2ang_0`, `wheel3ang_0`: 

    encode the deviation of the zero point of the telescope's motors:
    - _wheel1_ correspond to the boresight motor, and the zero point angle cause a rotation around the z-axis;
    - _wheel2_ correspond to the altitude motor, and the zero point angle cause a rotation around the y-axis;
    - _wheel3_ correspond to the ground motor, and the zero point angle cause a rotation around the z-axis.
    
2. `forkang`: 

    encodes the deviation of orthogonality between the H-AXIS and the V-AXIS;
    as reported in the fig. 2 this angles cause a rotation of the system around the x-axis.

3. `omegaVAXang`, `zVAXang`: 

    encode the deviation of the V-AXIS from the local vertical; 
    zVax is the angle between V-AXIS and the local vertical, 
    while omegaVAX is the azimuth of the ascending node (see fig. 3)

These angles are defined using a [`TelescopeAngles`](@ref) struct.

```@raw html
<figure>
    <img src="../assets/prm/fork.png" width="50%"/>
    <figcaption>Fig.2 Fork angle</figcaption>
</figure>
<figure>
    <img src="../assets/prm/wobble.png" width="50%"/>
    <figcaption>Fig.3 Wobble angles</figcaption>
</figure>
```
## Camera angles

To describe the orientation of a specific detector into the telescope reference frame we need to use three
Tait-Bryan angles:

- `panang`: representing a rotation around the x-axis

- `tiltang`: representing a rotation around the y-axis

- `rollang`: representing a rotation around the z-axis

With this convention a roll rotates the image seen by the camera around its center, a small pan (or tilt) shifts 
the image along the camera X (or Y) axis.

These angles are defined using a [`CameraAngles`](@ref) struct.

## Pointing Reconstruction Method

PRM consist of calculate a chain of rotations to project the direction of sight of a generic camera
into the Topocentric Horizontal Reference Frame (the ground r.f. of the telescope 
see [`telescopetoground`](@ref)).

For clarity, we can split the rotations in three steps ($R_i$ is a quaternion or a rotation matrix representing 
the rotation around the i axis):

$\mathbf{R}^{(\mathrm{tel})} = R_x(\mathrm{panang})R_y(\mathrm{tiltang})R_z(\mathrm{rollang})$
$\mathbf{R}^{(\mathrm{V-AXIS})} = R_z(\mathrm{wheel3ang-wheel3ang_0})R_x(\mathrm{forkang})R_y(\mathrm{wheel2ang-wheel2ang_0})$
$\mathbf{R}^{(\mathrm{geo})} = R_z(\mathrm{\omega_{VAX}})R_x(\mathrm{z_{VAX}})$

Taking as a reference Fig.1: the first step project the coordinates from the camera reference frame to the telescope r.f.,
then the second one from the telescope r.f. to the V-AXIS r.f. and at the last one from the V-AXIS r.f to the ground r.f.;
the final rotation operator describing the projection of the coordinate axis of the camera reference frame into the local 
topocentric reference frame is: $\mathbf{R}^{(\mathrm{geo})} * \mathbf{R}^{(\mathrm{V-AXIS})} * \mathbf{R}^{(\mathrm{tel})}$.

The $\mathbf{R}^{(\mathrm{tel})}$ operator is the output of [`camtotelescope`](@ref) while the other two operators are
computed by [`telescopetoground`](@ref).

## Example

Let's see a simple example, we'll now simulate the region of the sky observed by the telescope with and without the configuration
angles (we will use very large values i.e. not a very realistic case). Just like for [this example](@ref scanning_example) let's load all the
needed packages:

```@example prm
using Plots
using Healpix
using Stripeline
```

We can now define the control angles (representing the motor position in function of time) and the configuration angles ([`TelescopeAngles`](@ref) and 
[`CameraAngles`](@ref)) representing the non idealities of the system:

```@example prm
# Control angles
telescope_motors(time_s) = (0.0, deg2rad(20.0), timetorotang(time_s, 1))

# Configuration angles
telescope_ang = TelescopeAngles(
    forkang_rad = deg2rad(13.0),
    zVAXang_rad = deg2rad(10.0),
    omegaVAXang_rad = deg2rad(15.0)
)

detector_ang = CameraAngles()
nothing; #hide
```
For this example we will use the default values for the camera angles, meaning that the detector point towards
the boresight direction [0.0,0.0,1.0] perpendicular to the focal plane.
Now we can define a function that call [`genpointings`](@ref) with and without the
non idealities, iterate over matrix containing the directions and set a specific value for each pixel in a Healpix map:

```@example prm
function project_to_map(time_range, map, detector_ang, telescope_ang)
    # Call genpointings for the ideal case
    dirs, _ = genpointings(
        telescope_motors,
        detector_ang,
        time_range,
    )

    # Call genpointings for the non ideal case
    dirs_nonideal, _ = genpointings(
        telescope_motors,
        detector_ang,
        time_range,
        telescope_ang = telescope_ang
    )

    # For each sample, set the corresponding pixel in the sky map to:
    # 1 for the ideal case
    # 2 for the non ideal case
    for idx in 1:length(time_range)
        colat, long = dirs[idx, :]
        pixel_index = ang2pix(map, colat, long)
        map[pixel_index] = 1
        # Set the non ideal directions pixel
        colat, long = dirs_nonideal[idx, :]
        pixel_index = ang2pix(map, colat, long)
        map[pixel_index] = 2
    end
end
nothing; #hide
```

Finally, we can create the map calling `project_to_map` and plotting the result:

```@example prm
map = HealpixMap{Float64, RingOrder}(128)
sampling_time_s = 0.05
project_to_map(0.0:sampling_time_s:60.0, map, detector_ang, telescope_ang)
plot(map, orthographic)
savefig("oneminutemap_prm.svg"); nothing # hide
```

![](oneminutemap_prm.svg)

Where the pixels set to 2 are the non ideal case taking into account of
the configuration angles (in this example only the _forkang_ and the
_wobble angles_), while the pixels set to 1 are the ideal case alredy
discussed [here](@ref scanning_example).

## Reference Documentation

For a complete list of function used to reconstruct the scanning
direction, see [this](@ref scanning_docs).  

```@docs
TelescopeAngles
CameraAngles
camtotelescope
```