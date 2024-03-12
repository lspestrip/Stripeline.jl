export StripTod, StokesTod, split_into_n, allocate_tod
export stokes, ntimelines
export pwrq1, pwrq2, pwru1, pwru2
export demq1, demq2, demu1, demu2
export stokesi, stokesq, stokesu
export PWR_Q1_RANK, PWR_Q2_RANK, PWR_U1_RANK, PWR_U2_RANK
export DEM_Q1_RANK, DEM_Q2_RANK, DEM_U1_RANK, DEM_U2_RANK
export STOKES_I_RANK, STOKES_Q_RANK, STOKES_U_RANK
export fillnoise!

import Printf: @printf
import Random: AbstractRNG
import FFTW: fft, ifft, fftfreq
import LinearAlgebra: rmul!

using AstroTime
using RandomNumbers.PCG

"""
    split_into_n(length, num_of_segments)

Given the `length` of the array, the function returns a split of the
elements in `num_of_segments` subsets whose length is roughly the
same.

The result is an array containing the number of elements which must go
in each section; all the elements in the array have the property of
being similar and their sum being equal to `length`.

# Examples
```julia-repl
julia> split = split_into_n(20, 3)
3-element Array{Int64,1}:
6
7
7

julia> sum(split)
20
```
"""
function split_into_n(length, num_of_segments)
    @assert num_of_segments > 0
    @assert length >= num_of_segments
    start_pos = zeros(Int, num_of_segments + 1)

    for i = 1:num_of_segments+1
        start_pos[i] = floor(((i - 1) * length / num_of_segments))
    end

    return start_pos[2:end] - start_pos[1:end-1]
end


abstract type AbstractTod end

ntimelines(t::Type{AbstractTod}) = 0  # This makes allocate_tod fail

@doc raw"""
A matrix of time-ordered PWR/DEM data for a subset of polarimeters

The fields of this structure are the following:

- `time_range`: a range representing the time of each sample in the
  TOD.

- `samples`: a N×8×m matrix, where N is the number of samples in the
  timeline, and m is the number of polarimeters. The middle dimension
  spans the three Stokes parameters `I`, `Q`, and `U` (in this order).
  The value `N` is equal to `length(time_range)`

- `rng`: a pseudo-random number generator to be used for this TOD. It
  is guaranteed that the generator is uncorrelated with any other
  generator used for other `StripTod` objects created by the same call
  to `allocate_tod`.

"""
struct StripTod{T<:Number,S,R<:AbstractRNG} <: AbstractTod
    polarimeters::Any
    time_range::AbstractRange{S}
    samples::Array{T,3}
    rng::R
end

ntimelines(t::Type{StripTod}) = 8 # PWR Q1/Q2/U1/U2 + DEM Q1/Q2/U1/U2

function Base.show(io::IO, tod::StripTod{T,S}) where {T,S}
    print(
        io,
        string(
            "TOD($(length(tod.polarimeters)) polarimeters, ",
            "$(length(tod.time_range)) rows/polarimeter, ",
            "time range from $(minimum(tod.time_range)) to ",
            "$(maximum(tod.time_range)))",
        ),
    )
end


@doc raw"""
A matrix of time-ordered PWR/DEM data for a subset of polarimeters

The fields of this structure are the following:

- `time_range`: a range representing the time of each sample in the
  TOD.

- `samples`: a N×8×m matrix, where N is the number of samples in the
  timeline, and m is the number of polarimeters. The middle dimension
  spans the three Stokes parameters `I`, `Q`, and `U` (in this order).
  The value `N` is equal to `length(time_range)`

- `rng`: a pseudo-random number generator to be used for this TOD. It
  is guaranteed that the generator is uncorrelated with any other
  generator used for other `StripTod` objects created by the same call
  to `allocate_tod`.

"""
struct StokesTod{T<:Number,S,R<:AbstractRNG} <: AbstractTod
    polarimeters::Any
    time_range::AbstractRange{S}
    samples::Array{T,3}
    rng::R
end

ntimelines(t::Type{StokesTod}) = 3  # I, Q, U

function Base.show(io::IO, tod::StokesTod{T,S}) where {T,S}
    print(
        io,
        string(
            "StokesTOD($(length(tod.polarimeters)) polarimeters, ",
            "$(length(tod.time_range)) rows/polarimeter, ",
            "time range from $(minimum(tod.time_range)) to ",
            "$(maximum(tod.time_range)))",
        ),
    )
end


const PWR_Q1_RANK = 1
const PWR_Q2_RANK = 2
const PWR_U1_RANK = 3
const PWR_U2_RANK = 4
const DEM_Q1_RANK = 5
const DEM_Q2_RANK = 6
const DEM_U1_RANK = 7
const DEM_U2_RANK = 8
const STOKES_I_RANK = 1
const STOKES_Q_RANK = 2
const STOKES_U_RANK = 3

for (symb, ch_name, value) in [
    (:PWR_Q1_RANK, "PWR_Q1", PWR_Q1_RANK),
    (:PWR_Q2_RANK, "PWR_Q2", PWR_Q2_RANK),
    (:PWR_U1_RANK, "PWR_U1", PWR_U1_RANK),
    (:PWR_U2_RANK, "PWR_U2", PWR_U2_RANK),
    (:DEM_Q1_RANK, "DEM_Q1", DEM_Q1_RANK),
    (:DEM_Q2_RANK, "DEM_Q2", DEM_Q2_RANK),
    (:DEM_U1_RANK, "DEM_U1", DEM_U1_RANK),
    (:DEM_U2_RANK, "DEM_U2", DEM_U2_RANK),
    (:STOKES_I_RANK, "STOKES_I", STOKES_I_RANK),
    (:STOKES_Q_RANK, "STOKES_Q", STOKES_Q_RANK),
    (:STOKES_U_RANK, "STOKES_U", STOKES_U_RANK),
]
    @eval begin
        @doc """
    $($ch_name)_RANK = $($value)

The rank of `$($ch_name)` in the `samples` field of a [`StripTod`](@ref) variable.

"""
        $symb
    end
end


for (index, ch_name, fn_name) in [
    (PWR_Q1_RANK, "PWR_Q1", :pwrq1),
    (PWR_Q2_RANK, "PWR_Q2", :pwrq2),
    (PWR_U1_RANK, "PWR_U1", :pwru1),
    (PWR_U2_RANK, "PWR_U2", :pwru2),
    (DEM_Q1_RANK, "DEM_Q1", :demq1),
    (DEM_Q2_RANK, "DEM_Q2", :demq2),
    (DEM_U1_RANK, "DEM_U1", :demu1),
    (DEM_U2_RANK, "DEM_U2", :demu2),
    (STOKES_I_RANK, "STOKES_I", :stokesi),
    (STOKES_Q_RANK, "STOKES_Q", :stokesq),
    (STOKES_U_RANK, "STOKES_U", :stokesu),
]
    @eval begin
        function ($fn_name)(t::StripTod, pol::T) where {T<:Integer}
            @view t.samples[:, $index, pol]
        end

        @doc """
    $($fn_name)(tod::Stripeline.StripTod, pol::T) where T <: Integer

Return a view of the timeline for channel `$($ch_name)` in the
time-ordered data for the polarimeter indexed by `pol` in `tod`. It is
a shortcut for

    tod.samples[:, $($ch_name)_RANK, pol]

"""
        $fn_name
    end
end


@doc raw"""
    allocate_tod(todtype, t, time_range, polarimeters; mpi_rank=0, mpi_size=1, zero_tod=true)

Allocate memory for a time-ordered data matrix appropriate for the set
of polarimeters in the variable `polarimeters` (an array or range,
only its length is used and it is copied as-is in the returned
[`StripTod`](@ref) value). The type of each sample is set to `t`. The
variable `time_range` must be a range of floating-point values or
epochs (as defined by the
[AstroTime](https://juliaastro.github.io/AstroTime.jl/stable/)
package) which represent the full set of samples acquired during an
observation.

The optional parameters `mpi_rank` and `mpi_size` can be used in MPI
scripts to allocate only the subset of samples that must be processed
by the current MPI process.

If `zero_tod` is true, each sample in the TOD is set to zero;
otherwise, its value is undefined. Using `zero_tod = false` is
slightly faster, but it should be used only if you are going to set by
hand every sample in the TOD once `allocate_tod` returns.

Return a [`StripTod`](@ref) value.

### Example ###

```jldoctest; setup = :(using Stripeline)
julia> tod = allocate_tod(StripTod, Float32, 0.0:0.01:100.0, 1:5)
TOD(5 polarimeters, 10001 rows/polarimeter, time range from 0.0 to 100.0)
```

"""
function allocate_tod(
    todtype::Type{Tod},
    ::Type{T},
    time_range::R,
    polarimeters;
    mpi_rank = 0,
    mpi_size = 1,
    zero_tod = true,
    rng_seed = 12345,
) where {Tod<:AbstractTod,S,V,R<:AbstractRange,T<:Real}

    rank_idx = mpi_rank + 1 # Same as rank, but it starts from 1

    nsamples_arr = split_into_n(length(time_range), mpi_size)
    nsamples = nsamples_arr[rank_idx]

    first_idx = sum(nsamples_arr[1:(rank_idx-1)]) + 1
    last_idx = first_idx + nsamples - 1

    # Determine the time range for the current rank
    cur_range = time_range[first_idx]:step(time_range):time_range[last_idx]

    # Size of the field `samples`
    matrix_size = (length(cur_range), ntimelines(todtype), length(polarimeters))

    master_rng = PCGStateOneseq(UInt64, rng_seed)
    # Move the RNG seed to the position of this MPI rank
    advance!(master_rng, rank_idx)

    # Get the pseudorandom seed for this TOD
    cur_rng_seed = rand(master_rng, UInt64)
    rng = PCGStateSetseq(UInt64, PCG_XSH_RR, (cur_rng_seed, UInt64(rank_idx)))

    todtype(
        polarimeters,
        cur_range,
        zero_tod ? zeros(T, matrix_size) : Array{T}(undef, matrix_size),
        rng,
    )
end


@doc raw"""
    stokes(tod::StripTod{T, S, R})

Return a [`StokesTod`](@ref) variable containing the estimated
``I``/``Q``/``U`` Stokes parameters associated with the parameter
`tod`, which is a time-ordered data variable containing the eight
outputs `PWR_Q1`, `PWR_Q2`, `PWR_U1`, `PWR_U2`, `DEM_Q1`, `DEM_Q2`,
`DEM_U1`, `DEM_U2` of a set of Strip polarimeters.

"""
function stokes(tod::StripTod{T,S,R}) where {T,S,R}
    nelems, _, npols = size(tod.samples)
    samples = Array{T,3}(undef, (nelems, ntimelines(StokesTod), npols))

    for polidx = 1:npols
        samples[:, STOKES_I_RANK, polidx] =
            (
                pwrq1(tod, polidx) +
                pwrq2(tod, polidx) +
                pwru1(tod, polidx) +
                pwru2(tod, polidx)
            ) / 4
        samples[:, STOKES_Q_RANK, polidx] = (demq1(tod, polidx) + demq2(tod, polidx)) / 2
        samples[:, STOKES_U_RANK, polidx] = (demu1(tod, polidx) + demu2(tod, polidx)) / 2
    end

    return StokesTod(tod.polarimeters, tod.time_range, samples, tod.rng)
end


@doc raw"""

    fillnoise!(covfn_k2, tod::StripTod)

Overwrite all the samples in the TOD with white noise so that the
covariance matrix of the PWR and DEM signals for each polarimeter
is determined by the result of a call to `covfn_k2` like the following:

    (pwr_cov, dem_cov) = covfn_k2(polarimeter)

where `polarimeter` is one of the values in `tod.polarimeters`.
The two variables `pwr_cov` and `dem_cov` must be 4×4 covariance
matrices that express the covariance in K² (**not** in ADU!).

A typical call to `fillnoise!` will use the instrument database to
retrieve the covariance matrices; for instance:

    db = Sl.InstrumentDB()
    tod = Sl.allocate_tod(Float32, 0.0:0.1:10000.0, ["I3"])
    Tmp.fillnoise!(tod) do polname
        spec = Sl.spectrum(db, polname)
        (spec.pwr_cov_matrix_k2, spec.dem_cov_matrix_k2)
    end

No 1/f noise is simulated in this function.

"""
function fillnoise!(covfn_k2, tod::StripTod{T,S,R}) where {T,S,R}
    # Generate white noise
    for i in eachindex(tod.samples)
        tod.samples[i] = randn(tod.rng)
    end

    # Correlate the noise and scale it
    for (polidx, pol) in enumerate(tod.polarimeters)
        pwr_cov, dem_cov = covfn_k2(pol)

        pwr_chol = cholesky(Matrix(pwr_cov))
        pwr_samples = @view tod.samples[:, PWR_Q1_RANK:PWR_U2_RANK, polidx]
        rmul!(pwr_samples, pwr_chol.U)

        dem_chol = cholesky(Matrix(dem_cov))
        dem_samples = @view tod.samples[:, DEM_Q1_RANK:DEM_U2_RANK, polidx]
        rmul!(dem_samples, dem_chol.U)
    end
end
