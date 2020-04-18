```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Simulating the scanning strategy

Stripeline includes a set of tools to simulate the observation of the
sky by the telescope. The way a telescope observes the sky is called
*scanning strategy*, and it is obviously one of the most basic and
important task in any simulation of data taking.

## Strip's location and movements

The Strip telescope, unlike most of the optical telescopes, performs a
regular, uninterrupted movement using the so-called *azimuth motor*, a
wheel that makes the whole structure spin around its gravity axis. The
height of the telescope can be varied by another motor (the *altitude
motor*), but we foresee that it will be always kept at the same angle
during spinning. The following animation shows how the two motor
operates: first, the altitude motor moves by roughly 20°, and then the
azimuth motor starts spinning.

```@raw html
<iframe
    src="https://www.youtube.com/embed/01P_5aCmST0"
    frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
    allowfullscreen>
</iframe>
```

The red laser beam shows the so-called *boresight direction*, which
points toward the position being observed on the celestial sphere. The
motion of the telescope, coupled with the spinning of the Earth,
produces a characteristic pattern in the sequence of points on the
celestial sphere that are observed by Strip. To understand this
pattern, you must be aware that Strip will be deployed to Tenerife,
one of the Canarian Islands (Spain), a location in the Northern
Emisphere, at latitude ≈28°N:

```@raw html
<iframe
    width="425"
    height="350"
    frameborder="0"
    scrolling="no"
    marginheight="0"
    marginwidth="0"
    src="https://www.openstreetmap.org/export/embed.html?bbox=-16.525948047637943%2C28.29287470104711%2C-16.497623920440677%2C28.309538133052932&amp;layer=mapnik"
    style="border: 1px solid black">
</iframe>
<br/>
<small>
    <a href="https://www.openstreetmap.org/#map=16/28.3012/-16.5118">Enlarged map</a>
</small>
```

(The coordinates of the site, as well as its elevation, are available
in the constants [`TENERIFE_LATITUDE_DEG`](@ref),
[`TENERIFE_LONGITUDE_DEG`](@ref), and [`TENERIFE_HEIGHT_M`](@ref).)
Now, consider the animation of the spinning telescope as seen from
space:

```@raw html
<iframe
    src="https://www.youtube.com/embed/WKibmFWynbM"
    frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
    allowfullscreen>
</iframe>
```

As above, the red vector shows the boresight direction, while the
green vector shows the zenith at Tenerife (the direction along the
gravitational force vector, pointings towards the sky); the yellow
direction represents the direction of the beam polarizer, and it is a
parameter used to convert the [Stokes
parameters](https://en.wikipedia.org/wiki/Stokes_parameters) ``Q`` and
``U`` to celestial coordinates. Remembering that the red beam shows
where the telescope is looking, this scanning strategy observes a
circumference in the sky. However, because Earth rotates once every 24
hours (much faster in the animation than in reality!), the observation
pattern is a spiral that covers a «strip» on the sky, whose apparent
height is roughly equal to twice the angle of the boresight wheel.

## Simulating the motion of the telescope wheels

We'll now use the tools provided by Stripeline to simulate the region
of the sky that is observed by the telescope. You can type the
following commands in a Jupyter notebook, if you have installed the
awesome [IJulia](https://github.com/JuliaLang/IJulia.jl) package;
otherwise, you can write them directly on Julia's command line. Let's
load the packages we'll need to run the simulation:
[Plots](https://github.com/JuliaPlots/Plots.jl) enables the `plot`
function, [Healpix](https://github.com/ziotom78/Healpix.jl) implements
the Healpix tessellation algorithm on the sphere, and yours truly
Stripeline:

```@example scanningstrategy
using Plots
using Healpix
using Stripeline
```

Since Stripeline is a simulator, it has been designed to be extremely
versatile: the angles of the two wheels discussed above (the elevation
wheel and the azimuth wheel) can be specified as arbitrary functions
of time, expressed in seconds:

```@example scanningstrategy
telescope_motors(time_s) = (0.0, deg2rad(20.0), timetorotang(time_s, 1))
nothing; # hide
```

Here we use the handy function [`timetorotang`](@ref) to compute
``2\pi\nu t``, the angle as a function of time for a continuous
circular motion of 1 rotation per minute. This function must be passed
as the first parameter to the function [`genpointings`](@ref), which
computes the direction and the observing angle for a number of points
in time. The direction is encoded as a ``N\times 2`` matrix containing
the colatitude and the longitude in the two columns; we iterate on the
``N`` rows of this matrix and set each pixel in a Healpix map to this
value, taking advantage of Healpix.jl's function `ang2pix`:

```@example scanningstrategy
function project_to_map(time_range, map)
    # We discard the second return value of "genpointings" with "_":
    # it's the orientation, but for this simple example it is not necessary
    dirs, _ = genpointings(
        telescope_motors,  # Angle of the three wheels as a function of time
        Float64[0, 0, 1],  # Observing direction in the focal plane reference frame
        time_range,        # Array of time values
    )

    # For each sample, set the corresponding pixel in the sky map to 1
    for idx in 1:length(time_range)
        # Extract the colatitude and the longitude from "dirs"
        colat, long = dirs[idx, :]
        pixel_index = ang2pix(map, colat, long)
        map[pixel_index] = 1
    end
end
nothing; # hide
```

It's now a matter of creating a Healpix map, wrapping the code in a
call to `project_to_map`, and plotting the result:

```@example scanningstrategy
# We create a Healpix map that represents the whole sky sphere,
# tessellated up to some resolution NSIDE=128
map = Healpix.Map{Float64, RingOrder}(128)

# Set the sampling time, the number of seconds used to create one
# measurement. The real value is 0.01, but we use a larger
# sampling time to make the simulation faster
sampling_time_s = 0.05

# Run a simulation lasting one minute (which corresponds to
# exactly one rotation of the telescope around the azimuth axis)
project_to_map(0.0:sampling_time_s:60.0, map)

# Plot the map, using the orthographic projection
plot(map, orthographic)
savefig("oneminutemap.svg"); nothing # hide
```

![](oneminutemap.svg)

The result should match your expectations: considering the animation
of the spinning telescope as seen from space, the sequence of points
set to 1 in the map corresponds to the directions visited by the red
vector during *one* spin. Let's re-run the animation over a longer
time span, one hour; no need to re-allocate the map, as we'll simply
overwrite it:

```@example scanningstrategy
# Run a simulation lasting one hour (60 rotations
# of the telescope, while the Earth slowly spins)
project_to_map(0.0:sampling_time_s:3600.0, map)

plot(map, orthographic)
savefig("onehourmap.svg"); nothing # hide
```

![](onehourmap.svg)

Running the simulation over one hour shows the effect of the Earth's
rotation, which makes the circle move Eastward.

## A more complex example

In the examples above, we set to 1 the pixels in the map that were
observed by the Strip telescope at least one time. A much more
informative way of plotting these maps is to *count* the number of
times a pixel has been observed, as typically one wants to make
several observations and then average them together to produce one
estimate with its own error bar. The number stored in each pixel is
the so-called *hit count* for that pixel, and it is related to the
overall amount of time spent by the telescope observing that
direction.

Let's modify the function `project_to_map` so that it increments the
value of a pixel:

```@example scanningstrategy
function project_to_map(time_range, map)
    dirs, _ = genpointings(telescope_motors, Float64[0, 0, 1], time_range)

    for idx in 1:length(time_range)
        colat, long = dirs[idx, :]
        pixel_index = ang2pix(map, colat, long)

        # Increment the "hit count"
        map[pixel_index] += 1
    end
end
nothing; # hide
```

Now let's recreate the simulation above, using the `@animate` macro to
produce a movie of the hit count as time passes:

```@example scanningstrategy
# Reset the map, so that each pixel is set to zero
map[:] .= 0

# Create one frame per each minute of observation
anim = @animate for minute in 0:60
    # Start and end times of this minute; note that we drop
    # the last sample from "end_time_s", as it will be included
    # in the next iteration
    start_time_s = minute * 60
    end_time_s = (minute + 1) * 60 - sampling_time_s

    project_to_map(start_time_s:sampling_time_s:end_time_s, map)

    # The keyword "clim" fixes the lower and upper limits of the
    # color bar. Avoiding to do so results in the color scale flickering
    # between frames in the animation (try it!)
    plot(map, orthographic, clim=(0, 100))
end

# Save the result into an animated GIF file, with 10 frames per second.
# As one frame in our animation lasts one minute, this means that each
# second in the animation corresponds to 10 minutes
gif(anim, "scanning-animation.gif", fps = 10)
```

Note that most of the pixels are observed a few times, but those on
the uppermost and lowermost part of the strip have a much higher hit
count. We clipped the maximum value shown in the color bar to 100, but
we can make Julia compute the maximum value in the map:

```@example scanningstrategy
println("The maximum hit count in the map is ", maximum(map))
```

## A more complex example

So far, we have simulated the behavior of the Strip instrument in
quite ideal conditions: we started each simulation from time ``t =
0``, without specifying an absolute date.

The function [`genpointings`](@ref) accepts a starting time expressed
as a [Julian date](https://en.wikipedia.org/wiki/Julian_day); in this
case, it uses a slower algorithm to find the direction and orientation
of the telescope that considers several other effects:

1.  The position of the Earth with respect to the Sun;
2.  The [precession of the Earth's
    axis](https://en.wikipedia.org/wiki/Axial_precession);
3.  The
    [nutation](https://en.wikipedia.org/wiki/Astronomical_nutation) of
    the Earth's axis due to other bodies in the Solar System;
4.  [Stellar aberration](https://en.wikipedia.org/wiki/Aberration_(astronomy));
5.  [Atmospheric
    refraction](https://en.wikipedia.org/wiki/Atmospheric_refraction),
    although this computation is valid only for visible wavelengths
    and should therefore not be used when simulating observations done
    by Strip detectors (which operates at microwave lengths).

To specify times, you can use the functions `jdcnv` and `daycnv` from
[AstroLib](https://github.com/JuliaAstro/AstroLib.jl):

```@example scanningstrategy
using AstroLib, Dates

# Assume that the simulation starts on January, 1st 2022, 15:00:00
start_day = DateTime(2022, 1, 1, 15, 0, 0))

println("The simulation starts from JD $start_day")

dirs, orientations = genpointings(
    telescope_motors,
    Float64[0, 0, 1],
    0:sampling_time_s:60.0,
    start_day,   # We specify here the JD when the simulation starts
)
```

Passing `start_day` will make [`genpointings`](@ref) use a much slower
algorithm to compute the pointings; you should use this syntax only if
your simulation really needs the increased precision. Typical examples
where you really want to do this are the following:

1.  You want to estimate when and how long an object in the sky (e.g.,
    the [Crab Nebula](https://en.wikipedia.org/wiki/Crab_Nebula),
    [Jupiter](https://en.wikipedia.org/wiki/Jupiter)) will be visible
    by Strip;
2.  A variation of the previous point is to compute when bright
    objects (Sun, Moon, etc.) might be dangerously close to the
    boresight direction of the telescope, in order to decide when
    the telescope will need to be shut down to prevent overheating or
    saturations in the detectors.
3.  You want to produce a sky map to be compared with those of other
    experiments running at the same time as Strip;
4.  You want to assess how much effects like precession, nutation, and
    aberration affect the measurements.

## Reference documentation

```@docs
TENERIFE_LATITUDE_DEG
TENERIFE_LONGITUDE_DEG
TENERIFE_HEIGHT_M
telescopetoground
groundtoearth
genpointings
timetorotang
northdir
eastdir
polarizationangle
```
