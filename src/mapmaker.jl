export TodNoiseProperties, build_noise_properties, tod2map_mpi, baseline2map_mpi, destripe, baselines_covmat

using LinearAlgebra

try
    import MPI
catch
end

@doc raw"""
This structure holds a number of parameters relative to the noise level and
the 1/f baselines measured or simulated for a certain polarimeter.

Field              | Type           | Meaning
:----------------- |:-------------- |:----------------------------------------------
`pol_number`       | Int            | ID number of the polarimeter
`sigma`            | Float          | sigma of white noise
`num_of_samples`   | Int            | total number of samples
`baseline_lengths` | Array{Int}     | Array containing the length of each baseline

The structure provides a constructor that accepts the following
three keyword parameters:

- `pol`: number of the polarimeter
- `rms`: the sigma level
- `baseline_lengths`: array containing the elements of each baseline

# Example

The followng code creates a `TodNoiseProperties` object referring to polarimeter #1
(int this case, STRIP31) and specifying 0.002463 as the noise level. The data sample
contains 500 samples, and it is split into 5 baselines of equal length.

```julia
prop = TodNoiseProperties(
    pol = 1,
    rms = 0.002463, 
    baselines = Int[100, 100, 100, 100, 100, 100],
)
```

"""
struct TodNoiseProperties
    polarimeter::Int
    sigma::Float64
    number_of_samples::Int
    baseline_lengths::Array{Int}

    TodNoiseProperties(;
        pol::Int,
        rms::Float64,
        baselines::Array{Int}) = new(pol, rms, sum(baselines), baselines)
end


@doc raw"""
    function build_noise_properties(detector_list, rms_list, num_of_baselines, num_of_samples) -> data_properties

This function builds a list of `TodNoiseProperties` object. These are
needed for desriping and calculation of `baselines_covmat`.

It requires in input 4 arrays containing:
- the ID number of the polarimeters that the current rank will simulate
- the white noise sigma for each polarimeter that the current rank will simulate
- the number of 1/f baselines for each polarimeter that the current rank will simulate
- the total number of samples for each polarimeter that the current rank will simulate

N.B. the ID numbers, the number of 1/f baselines and the total number of samples can be obtained using the function `get_chunk_properties`.

It returns an array of `TodNoiseProperties`, of length equal to the number of polarimeters simulated by current rank.

"""
function build_noise_properties(detector_list, rms_list, num_of_baselines, num_of_samples)
    @assert length(detector_list) == length(rms_list)
    @assert length(detector_list) == length(num_of_baselines)
    @assert length(detector_list) == length(num_of_samples)

    data_properties  = Array{TodNoiseProperties}(undef, length(detector_list))
    for i in 1:length(detector_list)
        # We use baselines of equal lengths
        baseline_lengths  = repeat([Int(num_of_samples[i] / num_of_baselines[i])], num_of_baselines[i])
        data_properties[i] = TodNoiseProperties(
            pol = detector_list[i],
            rms = rms_list[detector_list[i]],
            baselines = baseline_lengths,
        )
    end
    data_properties
end


@doc raw"""
    tod2map(pix_idx, tod, num_of_pixels; comm = nothing) -> binned_map
    tod2map(pix_idx, tod, num_of_pixels, data_properties; comm = nothing) -> binned_map

This function creates a binned map from the time-ordered data kept in
the array `tod`, assuming that each sample is observing the pixel
whose index is in `pix_idx`. The parameter `num_of_pixels` contains
the number of pixels in the Healpix map, and it is used as an upper
bound for the values in `pix_idx`. The parameter `comm` must be a MPI
communicator, or `nothing` if you are not using MPI.

If `comm` is not `nothing`, the function is parallelized using MPI, and
each process computes a map from its available data.  All partial maps
are then combined together with `MPI.Allreduce`. The function returns
an array containing the binned map.

If Array of structures `TodNoiseProperties`is passed to the function, the
output binned map will be a weighted binned map. Each sample is weighted
according to the inverse white noise variance sigma^2 of the corrispondingù
polarimeter. In this way, the less noisy polarimeters will count more in the
estimation of the map.

# Requirements
- The length of the arrays `pix_idx` and `tod` must be the same

"""
function tod2map_mpi(pix_idx, tod, num_of_pixels; comm = nothing, unseen = NaN)
    @assert length(pix_idx) == length(tod)

    T = eltype(tod)
    N = eltype(pix_idx)

    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(N, num_of_pixels)
    binned_map = zeros(T, num_of_pixels)
    hits = zeros(N, num_of_pixels)

    @inbounds for i in eachindex(pix_idx)
        partial_map[pix_idx[i]] += tod[i]
        partial_hits[pix_idx[i]] += 1
    end

    if comm != nothing
        binned_map = MPI.Allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.Allreduce(partial_hits, MPI.SUM, comm)
    else
        binned_map .= partial_map
        hits .= partial_hits
    end

    @inbounds for i in eachindex(binned_map)
        if (hits[i] > 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end

    binned_map
end

function tod2map_mpi(pix_idx, tod, num_of_pixels, data_properties; comm = nothing, unseen = NaN)
    @assert length(pix_idx) == length(tod)

    T = eltype(tod)

    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(T, num_of_pixels)
    binned_map = zeros(T, num_of_pixels)
    hits = zeros(T, num_of_pixels)

    start_idx = 1

    for j in eachindex(data_properties)                 #loop on detectors
        end_idx = start_idx + data_properties[j].number_of_samples - 1
        for i in start_idx:end_idx             #loop on samples
            partial_map[pix_idx[i]] += tod[i] * 1 / (data_properties[j].sigma)^2
            partial_hits[pix_idx[i]] += 1 / (data_properties[j].sigma)^2
        end
        start_idx += data_properties[j].number_of_samples
    end

    if comm != nothing
        binned_map = MPI.Allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.Allreduce(partial_hits, MPI.SUM, comm)
    else
        binned_map .= partial_map
        hits .= partial_hits
    end

    @inbounds for i in eachindex(binned_map)
        if (hits[i] > 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end

    binned_map
end

@doc raw"""
    baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels; comm = nothing)-> noise_map
    baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties; comm = nothing) -> noise_map

This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with
MPI.Allreduce.
The function returns an array containing the binned map.

If Array of structures `TodNoiseProperties`is passed to the function (instead of `baseline_lengths`) the output binned map will be a weighted binned map.
Each sample will be weighted according to the inverse white noise variance sigma^2 of the corrisponding polarimeter.
In this way, the less noisy polarimeters will count more in the estimation of the map.

# Requirements
- The length of `baselines` and `baseline_lengths` must be the same;
- The value `sum(baseline_lengths)` must be the same as the length of `pix_idx`.
"""
baseline2map_mpi


function baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels; comm = nothing,
                          unseen = NaN)

    T = eltype(baselines)
    N = eltype(pix_idx)

    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(N, num_of_pixels)
    noise_map = zeros(T, num_of_pixels)
    hits = zeros(N, num_of_pixels)

    startidx = 1

    for i in eachindex(baseline_lengths)
        endidx = baseline_lengths[i] + startidx - 1

        for j in startidx:endidx
            partial_map[pix_idx[j]] += baselines[i]
            partial_hits[pix_idx[j]] += 1
        end
        startidx += baseline_lengths[i]
    end

    if comm != nothing
        noise_map = MPI.Allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.Allreduce(partial_hits, MPI.SUM, comm)
    else

        noise_map .= partial_map
        hits .= partial_hits
    end

    @inbounds for i in eachindex(noise_map)
        if (hits[i] > 0)
            noise_map[i] /= hits[i]
        else
            noise_map[i] = unseen
        end
    end

    noise_map
end

function baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties, num_of_baselines; comm = nothing, unseen = NaN)

    T = eltype(baselines)

    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(T, num_of_pixels)
    noise_map = zeros(T, num_of_pixels)
    hits = zeros(T, num_of_pixels)

    startidx = 1
    baseline_idx = 1

    for det_idx in eachindex(data_properties)  #loop on detectors
        for sample_idx in eachindex(data_properties[det_idx].baseline_lengths)
            endidx = data_properties[det_idx].baseline_lengths[sample_idx] + startidx - 1

            @inbounds for j in startidx:endidx
                partial_map[pix_idx[j]] += baselines[baseline_idx] / (data_properties[det_idx].sigma)^2
                partial_hits[pix_idx[j]] += 1 / (data_properties[det_idx].sigma)^2
            end
            startidx += data_properties[det_idx].baseline_lengths[sample_idx]
            baseline_idx += 1
        end
    end

    if comm != nothing
        noise_map = MPI.Allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.Allreduce(partial_hits, MPI.SUM, comm)
    else

        noise_map .= partial_map
        hits .= partial_hits
    end

    @inbounds for i in eachindex(noise_map)
        if (hits[i] > 0)
            noise_map[i] /= hits[i]
        else
            noise_map[i] = unseen
        end
    end

    noise_map
end


function applyz_and_sum(pix_idx, tod, num_of_pixels, data_properties, num_of_baselines; comm = nothing, unseen = NaN)

    @assert length(tod) == length(pix_idx)

    baselines_sum = zeros(eltype(tod), num_of_baselines)

    binned_map = tod2map_mpi(pix_idx,
        tod,
        num_of_pixels,
        data_properties,
        comm = comm,
        unseen = unseen)

    baseline_idx = 1
    startidx = 1
    for det_idx in eachindex(data_properties)  #loop on detectors
        for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
            endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

            for j in startidx:endidx
                baselines_sum[baseline_idx] += (tod[j] - binned_map[pix_idx[j]]) * 1 / (data_properties[det_idx].sigma)^2
            end

            startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
            baseline_idx += 1
        end
    end

    baselines_sum
end


function applya(baselines, pix_idx, num_of_baselines, num_of_pixels, data_properties; comm = nothing, unseen = NaN)
    @assert length(baselines) == num_of_baselines

    baselines_sum = zeros(eltype(baselines), num_of_baselines)
    total_sum = zero(eltype(baselines))

    binned_map = baseline2map_mpi(pix_idx,
        baselines,
        num_of_pixels,
        data_properties,
        num_of_baselines,
        comm = comm,
        unseen = unseen)

    startidx = 1
    baseline_idx = 1

    for det_idx in eachindex(data_properties)
        for baseline_idx in eachindex(data_properties[det_idx].baseline_lengths)
            endidx = data_properties[det_idx].baseline_lengths[baseline_idx] + startidx - 1

            for j in startidx:endidx
                baselines_sum[baseline_idx] += (baselines[baseline_idx] - binned_map[pix_idx[j]]) * 1 / (data_properties[det_idx].sigma)^2
            end

            startidx += data_properties[det_idx].baseline_lengths[baseline_idx]
            baseline_idx += 1
        end
    end

    #needed to assure that sum(baselines)==0

    if comm != nothing
        total_sum = MPI.Allreduce([sum(baselines)], MPI.SUM, comm)[1]
    else
        total_sum = sum(baselines)
    end

    baselines_sum .+= total_sum
end


function mpi_dot_prod(x, y; comm = nothing)

    @assert eltype(x) == eltype(y)
    result = zero(eltype(x))
    local_sum::eltype(x) = dot(x, y)

    if comm != nothing
        result = MPI.Allreduce([local_sum], MPI.SUM, comm)[1]
    else
        result = local_sum
    end

    return result
end

mutable struct DestripingResults
    threshold::Float64
    max_iter::Int
    convergence_param_list::Array{Float64}
    best_iteration::Int
    best_baselines::Array{Float64}
    best_sky_map::Any
    baseline_history::Any

    DestripingResults() = new(
        1e-9,
        1000,
        Float64[],
        -1,
        Float64[],
        nothing,
        nothing,
    )
end


function conj_grad(
    results::DestripingResults,
    baselines_sum,
    pix_idx,
    tod,
    num_of_pixels,
    data_properties,
    num_of_baselines,
    rank;
    comm = nothing,
    save_baseline_history = false,
    callback = nothing,
)

    T = eltype(tod)
    N = eltype(pix_idx)

    baselines = Array{T}(undef, num_of_baselines)
    r = Array{T}(undef, num_of_baselines)
    r_next = Array{T}(undef, num_of_baselines)
    p = Array{T}(undef, num_of_baselines)
    Ap = Array{T}(undef, num_of_baselines)

    baselines .= 0    #starting baselines
    convergence_parameter = zero(T)
    rdotr = zero(T)
    rdotr_next = zero(T)

    best_convergence_parameter = zero(T)
    best_baselines = zeros(T, num_of_baselines)
    results.best_iteration = 0

    r = baselines_sum - applya(baselines, pix_idx, num_of_baselines, num_of_pixels, data_properties, comm = comm)  #residual
    p .= r
    rdotr = mpi_dot_prod(r, r, comm = comm)
    best_convergence_parameter = sqrt(rdotr)

    results.baseline_history = save_baseline_history ? [] : nothing
    results.best_baselines = Array{Float64}(undef, num_of_baselines)
    results.best_baselines .= baselines

    rdotr == 0 && return results

    results.convergence_param_list = Float64[]
    iter_idx = 0
    while true
        Ap = applya(p, pix_idx,  num_of_baselines, num_of_pixels, data_properties, comm = comm)

        rdotr = mpi_dot_prod(r, r, comm = comm)
        pdotAp = mpi_dot_prod(p, Ap, comm = comm)

        alpha = rdotr / pdotAp
        @. baselines += alpha * p
        save_baseline_history && push!(results.baseline_history, copy(baselines))

        @. r_next = r - alpha * Ap
        rdotr_next = mpi_dot_prod(r_next, r_next, comm = comm)
        convergence_parameter = sqrt(rdotr_next)
        push!(results.convergence_param_list, convergence_parameter)

        if (convergence_parameter < best_convergence_parameter)
            best_convergence_parameter = convergence_parameter
            results.best_baselines .= baselines
            results.best_iteration = iter_idx
        end

        isnothing(callback) || callback(
            iter_idx = iter_idx,
            max_iter = results.max_iter,
            convergence_parameter = convergence_parameter,
            convergence_threshold = results.threshold,
        )
        
        ((convergence_parameter < results.threshold) || (iter_idx > results.max_iter)) && break

        beta = rdotr_next / rdotr

        @. p = r_next + beta * p
        r .= r_next
        iter_idx += 1
    end
end


function destriped_map(baselines, pix_idx, tod, data_properties, num_of_pixels, num_of_baselines; comm = nothing, unseen = NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, data_properties, comm = comm) - baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties, num_of_baselines, comm = comm)
end


@doc raw"""
    destripe(pix_idx, tod, num_of_pixels, data_properties, rank; comm = nothing, threshold = 1e-9, max_iter = 10000, save_baseline_history = false) -> DestripingResults

This MPI based function creates a map from a TOD and removes both 1/f
and white noise, using the destriping technique.

The parameters passed to the function have the following meaning:

- `pix_idx`: array containing the indices of the pixels visited by the instrument

- `tod`: the values measured by the polarimeters for each pixel
  (either I, Q, or U)

- `num_of_pixels`: the number of pixels in the map to be
   produced. This is used as an upper limit for the values in `pix_idx`

- `data_properties`: an array of structures `TodNoiseProperties`,
   holding information on each simulated polarimeter noise level,
   number and length of 1/f baselines and total number of samples.
   It can be obtained by using function `build_noise_properties`.

- `rank`: the rank of the current MPI process

- `comm`: the MPI communicator object.

The following arguments are optional:

- `threshold` is used by the conjugated-gradient algorithm. When the
   residual error of the iteration goes below this value, the iteration
   stops. The smaller the value, the more accurate the solution.

- `max_iter` is the maximum number of iterations to be executed in the
   conjugated-gradient algorithm. If the algorithm does not converge
   within this number of iterations, the process will quit without
   having reached the convergence threshold (see the `threshold`
   keyword above).

- If `save_baseline_history` is `true`, the return value will contain the
  sequence of baselines tested by the CG algorithm. Each MPI process will
  hold its own baselines.

The function returns a `DestripingResults` object containings the destriped
map, the sequence of baselines, and other information describing the
convergence of the CG algorithm.

# Remarks
- The length of the arrays `pix_idx` and `tod` must be the same;
- If you do not specify `comm`, no MPI will be used
"""
function destripe(
    pix_idx,
    tod,
    num_of_pixels,
    data_properties,
    rank;
    comm = nothing,
    threshold = 1e-9,
    max_iter = 1_000,
    save_baseline_history = false,
    unseen = NaN,
    callback = nothing,
)

    @assert length(pix_idx) == length(tod)

    num_of_baselines = 0
    for i in 1:length(data_properties)
        num_of_baselines += length(data_properties[i].baseline_lengths)
    end

    baselines_sum = applyz_and_sum(pix_idx,
        tod,
        num_of_pixels,
        data_properties,
        num_of_baselines,
        comm = comm,
        unseen = unseen,
    )

    results = DestripingResults()
    results.threshold = threshold
    results.max_iter = max_iter

    conj_grad(
        results,
        baselines_sum,
        pix_idx,
        tod,
        num_of_pixels,
        data_properties,
        num_of_baselines,
        rank,
        comm = comm,
        save_baseline_history = save_baseline_history,
        callback = callback,
    )

    # once we have an estimate of the baselines, we can build the destriped map
    results.best_sky_map = destriped_map(results.best_baselines,
        pix_idx,
        tod,
        data_properties,
        num_of_pixels,
        num_of_baselines,
        comm = comm,
        unseen = unseen,
    )

    results
end


@doc raw"""
    baselines_covmat(polarimeters, sigma_k, baseline_length_s, fsamp_hz, total_time) -> covariance_matrix

This function produces the covariance matrix of the 1/f baselines computed by the destriper.
As an approximation, we ignore nondiagonal terms (i.e. the correlations between different baselines).
The function output is therefore an array corresponding to the covariance matrix diagonal.
Another assumption is that the white noise variance stays constant over a given baseline.

The baseline error is computed according to equation 29 in https://arxiv.org/abs/0904.3623

The parameters passed to the function have the following meaning:
    -`polarimeters`: the array of polarimeters ID numbers
    -`sigma_k`: the array of corresponding white noise sigma (in K)
    -`baseline_length_s` = the length (in s) of each 1/f baseline
    -`fsamp_hz` = the sampling frequency (in Hz)
    -`total_time` = the duration (in s) of the observation

The function firstly computes an Array of structures `TodNoiseProperties`,
which matches the white noise variance sigma with the number of 1/f baselines and their lengths, for each polarimeter.
Then it computes the covariance matrix according to these matches.
"""
function baselines_covmat(polarimeters, sigma_k, baseline_length_s, fsamp_hz, total_time)

    ALL_data_properties = Array{TodNoiseProperties}(undef, length(polarimeters))
    baselines_per_pol =  Int64(total_time / baseline_length_s)
    for i in 1:length(polarimeters)
        baseline_lengths  = repeat([baseline_length_s * fsamp_hz], baselines_per_pol)
        ALL_data_properties[i]  = TodNoiseProperties(pol = polarimeters[i], 
            rms = sigma_k[i],
            baselines = baseline_lengths,
        )
    end

    covariance_matrix = []
    for i in 1:length(ALL_data_properties)
        partial_covmat = Array{Float64}(undef, length(ALL_data_properties[i].baseline_lengths))

        for j in 1:length(ALL_data_properties[i].baseline_lengths)
            sigma = ALL_data_properties[i].sigma
            this_baseline_length = ALL_data_properties[i].baseline_lengths[j]

            partial_covmat[j] = sigma^2 / this_baseline_length
        end
        covariance_matrix = append!(covariance_matrix, partial_covmat)
    end
    return covariance_matrix
end
