export StripTod, split_into_n, allocate_tod
export pwrq1, pwrq2, pwru1, pwru2
export demq1, demq2, demu1, demu2
export PWR_Q1_RANK, PWR_Q2_RANK, PWR_U1_RANK, PWR_U2_RANK
export DEM_Q1_RANK, DEM_Q2_RANK, DEM_U1_RANK, DEM_U2_RANK

import Printf: @printf
import Random: AbstractRNG

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
    start_pos = zeros(Int, num_of_segments+1)
    
    for i in 1:num_of_segments+1
        start_pos[i] =  floor(((i-1)*length/num_of_segments))
    end
   
    return start_pos[2:end]-start_pos[1:end-1]    
end


@doc raw"""
A matrix of time-ordered data for a subset of polarimeters

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
struct StripTod{T <: Real, S, R <: AbstractRNG}
    polarimeters::Any
    time_range::AbstractRange{S}
    samples::Array{T,3}
    rng::R
end

function Base.show(io::IO, tod::StripTod{T, S}) where {T, S}
    println(io, string(
        "TOD($(length(tod.polarimeters)) polarimeters, ",
        "$(length(tod.time_range)) samples/polarimeter, ",
        "time range from $(minimum(tod.time_range)) to ",
        "$(maximum(tod.time_range)))",
    ))
end


const PWR_Q1_RANK = 1
const PWR_Q2_RANK = 2
const PWR_U1_RANK = 3
const PWR_U2_RANK = 4
const DEM_Q1_RANK = 5
const DEM_Q2_RANK = 6
const DEM_U1_RANK = 7
const DEM_U2_RANK = 8

for (symb, ch_name, value) in [
    (:PWR_Q1_RANK, "PWR_Q1", PWR_Q1_RANK),
    (:PWR_Q2_RANK, "PWR_Q2", PWR_Q2_RANK),
    (:PWR_U1_RANK, "PWR_U1", PWR_U1_RANK),
    (:PWR_U2_RANK, "PWR_U2", PWR_U2_RANK),
    (:DEM_Q1_RANK, "DEM_Q1", DEM_Q1_RANK),
    (:DEM_Q2_RANK, "DEM_Q2", DEM_Q2_RANK),
    (:DEM_U1_RANK, "DEM_U1", DEM_U1_RANK),
    (:DEM_U2_RANK, "DEM_U2", DEM_U2_RANK),
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
]
    @eval begin
        function ($fn_name)(t::StripTod, pol::T) where T <: Integer
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
    allocate_tod(t, time_range, polarimeters; mpi_rank=0, mpi_size=1, zero_tod=true)

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
julia> tod = allocate_tod(Float32, 0.0:0.01:100.0, 1:5)
TOD(5 polarimeters, 10001 samples/polarimeter, time range from 0.0 to 100.0)
```

"""
function allocate_tod(
    t::Type{T},
    time_range::R,
    polarimeters;
    mpi_rank = 0,
    mpi_size = 1,
    zero_tod = true,
    rng_seed = 12345,
) where {S, V, R <: AbstractRange, T <: Real}

    rank_idx = mpi_rank + 1 # Same as rank, but it starts from 1
    
    nsamples_arr = split_into_n(length(time_range), mpi_size)
    nsamples = nsamples_arr[rank_idx]

    first_idx = sum(nsamples_arr[1:(rank_idx - 1)]) + 1
    last_idx = first_idx + nsamples - 1

    # Determine the time range for the current rank
    cur_range = time_range[first_idx]:step(time_range):time_range[last_idx]

    # Size of the field `samples`
    matrix_size = (length(cur_range), 8, length(polarimeters))

    master_rng = PCGStateOneseq(UInt64, rng_seed)
    # Move the RNG seed to the position of this MPI rank
    advance!(master_rng, rank_idx)

    # Get the pseudorandom seed for this TOD
    cur_rng_seed = rand(master_rng, UInt64)
    rng = PCGStateSetseq(
        UInt64,
        PCG_XSH_RR,
        (cur_rng_seed, UInt64(rank_idx)),
    )
    
    StripTod(
        polarimeters,
        cur_range,
        zero_tod ? zeros(T, matrix_size) : Array{T}(undef, matrix_size),
        rng,
    )
end
