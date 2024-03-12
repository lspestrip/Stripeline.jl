```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Time-ordered data (TOD)

Stripeline contains functions that create the timeline of scientific
data as produced by each polarimeter and distribute it across several
parallel processes, potentially running through MPI.

A timeline of scientific data is generically called a *TOD*, and the
polarimeters implemented in Strip produce a complex flow of data, as
they output *eight* series:

-   Four series are called *`PWR` series*: `PWR_Q1`,
    `PWR_Q2`, `PWR_U1`, and `PWR_U2`;
-   Four series are called *`DEM` series*: `DEM_Q1`,
    `DEM_Q2`, `DEM_U1`, and `DEM_U2`.
    
Each of these series is sampled at 50 Hz and contains scientific data
(in the preliminary unit tests done in 2017, the sampling frequency
was 25 Hz; keep this in mind if you happen to replicate some returns
of those ancient tests!); the raw instrument represents them as
integer numbers, but Stripeline usually assumes they have already been
calibrated to Kelvin and therefore uses floating-point numbers.

A TOD is represented by the structure [`StripTod`](@ref), which
contains the following fields:

-   `polarimeters`: a generic sequence of values, which represent each
    polarimeter stored in the TOD. Very common choices for this field
    are a range (for example, `1:49` represents all the Q-band
    polarimeters) or a list of strings, e.g., `["I0", "I3", "V2"]`.
-   `time_range`: a range representing the time of each sample in the
    TOD. Very common choices for this are floating-point ranges like
    `0.0:(1 / 50):100.0` (100 s sampled at 50 Hz) and ranges of
    astronomical times, created using the
    [AstroTime](https://juliaastro.github.io/AstroTime.jl/stable/)
    package.
-   `samples`: this is the matrix containing the actual scientific
    samples. It has a shape of ``N \times 8 \times m``, where ``N`` is
    the number of samples along the time axis and must be equal to
    `length(time_range)`, and ``m`` is the number of polarimeters,
    which must be equal to `length(polarimeters)`. The middle
    dimension is 8 because there are eight timelines produced by each
    polarimeter (see above).

Allocating the TOD is a matter of calling the function
[`allocate_tod`](@ref):

```@example tods
using Stripeline

# One minute of data sampled at 50 Hz, 10 polarimeters
tod = allocate_tod(StripTod, Float32, 0.0:(1.0 / 50.0):60.0, 1:10)
```

In a MPI environment, you can pass the keywords `mpi_rank` and
`mpi_size` to [`allocate_tod`](@ref) and the function will only
allocate the subset of the overall time range that is assigned to the
running MPI process.

Once a [`StripTod`](@ref) variable is allocated, accessing the samples
in the TOD can be done through the field `samples`, which is a matrix:

```@example tods
println("The shape of the TOD is ", size(tod.samples))
```

Accessing the right column can be cumbersome, as there are eight of
them (`PWR_Q1`, `PWR_Q2`, …), so Stripeline provides the constants
[`PWR_Q1_RANK`](@ref), [`PWR_Q2_RANK`](@ref), … for this:

```@example tods
# Print the minimum sample acquired by PWR_Q1 for polarimeter number #4
println(minimum(tod.samples[:, PWR_Q1_RANK, 4]))
```

The result is zero because [`allocate_tod`](@ref) sets all the samples
to zero by default. An alternative approach is to use one of the
functions [`pwrq1`](@ref), [`pwrq2`](@ref), [`pwru1`](@ref),
[`pwru2`](@ref), [`demq1`](@ref), [`demq2`](@ref), [`demu1`](@ref),
[`demu2`](@ref), which return a *view* of the original array and can
thus be used both for reading and writing:

```@example tods
# Set one element of PWR_Q1 for polarimeter number #4
pwrq1(tod, 4)[142] = 1.36
println(
    "The maximum element of polarimeter #4 is now ", 
    maximum(pwrq1(tod, 4)),
)
```

```@docs
StripTod
allocate_tod
PWR_Q1_RANK
PWR_Q2_RANK
PWR_U1_RANK
PWR_U2_RANK
DEM_Q1_RANK
DEM_Q2_RANK
DEM_U1_RANK
DEM_U2_RANK
pwrq1
pwrq2
pwru1
pwru2
demq1
demq2
demu1
demu2
```


# I/Q/U TODs

The scientific information kept in a [`StripTod`](@ref) variable is
encoded in eight timelines (PWR_Q1`, `PWR_Q2`, `PWR_U1`, `PWR_U2`,
`DEM_Q1`, `DEM_Q2`, `DEM_U1`, `DEM_U2`), but most of the information
in these timelines is redundant and is essentially
correlated/uncorrelated noise.

The eight timelines can be combined to provide a direct estimate of
the four Stokes parameters ``I``, ``Q``, and ``U`` in the reference
frame of the detector, using these formulae:

```math
\begin{aligned}
I &= \frac{\text{PWRQ1} + \text{PWRQ2} + \text{PWRU1} + \text{PWRU2}}4,\\
Q &= \frac{\text{PWRQ1} - \text{PWRQ2}}2,\\
U &= \frac{\text{PWRU1} - \text{PWRU2}}2.
\end{aligned}
```

The three ``I``/``Q``/``U`` timelines can be stored in a
[`StokesTod`](@ref) variable, which behaves exactly like a
[`StripTod`](@ref) but only contains three timeseries instead of
eight. To allocate it, you can use [`allocate_tod`](@ref) again,
passing a [`StokesTod`](@ref) parameter:

```@example tods
stokestod = allocate_tod(StokesTod, Float32, 0.0:(1.0 / 50.0):60.0, 1:10)
```

Or you can compute it starting from a [`StokesTod`](@ref) variable:

```@example tods
stokestod = stokes(tod)
```

The latter call works out-of-the-box in a MPI environment, of course.
You can access the ``I``/``Q``/``U`` timelines either using the
explicit indexes [`STOKES_I_RANK`](@ref), [`STOKES_Q_RANK`](@ref),
[`STOKES_U_RANK`](@ref), or the three functions [`stokesi`](@ref),
[`stokesq`](@ref), [`stokesu`](@ref), similarly to what was described
above for [`StripTod`](@ref).

```@docs
StokesTod
stokes
STOKES_I_RANK
STOKES_Q_RANK
STOKES_U_RANK
stokesi
stokesq
stokesu
```
