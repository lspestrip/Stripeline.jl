export tod2map_mpi, baseline2map_mpi, destripe

using LinearAlgebra

try
    import MPI
catch
end 

@doc raw"""
    tod2map(pix_idx, tod, num_of_pixels, comm) -> binned_map

This function creates a binned map from the time-ordered data kept in
the array `tod`, assuming that each sample is observing the pixel
whose index is in `pix_idx`. The parameter `num_of_pixels` contains
the number of pixels in the Healpix map, and it is used as an upper
bound for the values in `pix_idx`. The parameter `comm` must be a MPI
communicator, or `missing` if you are not using MPI.

This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with
MPI.allreduce.

The function returns an array containing the binned map.

# Requirements

- The length of the arrays `pix_idx` and `tod` must be the same;

"""
function tod2map_mpi(pix_idx, tod, num_of_pixels, comm=missing;
                     unseen=NaN)   
    T = eltype(tod)
    N = eltype(pix_idx)
    
    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(N, num_of_pixels)
    binned_map = zeros(T, num_of_pixels)
    hits = zeros(N, num_of_pixels )

    for i in eachindex(pix_idx)
        partial_map[pix_idx[i]] += tod[i]
        partial_hits[pix_idx[i]] += 1
    end

    if(!ismissing(comm))
        binned_map = MPI.allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.allreduce(partial_hits, MPI.SUM, comm)
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
    baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm; unseen=NaN)

This function creates a binned map from the sequence of baselines in
`baselines`. Each baseline covers a number of samples equal to the
corresponding element in `baseline_lengths`. The function assumes that
each sample is observing the pixel whose index is in `pix_idx`. The
parameter `num_of_pixels` contains the number of pixels in the Healpix
map, and it is used as an upper bound for the values in `pix_idx`. The
parameter `comm` must be a MPI communicator, or `missing` if you are
not using MPI.

This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with
MPI.allreduce.

The function returns an array containing the binned map.

# Requirements

- The length of `baselines` and `baseline_lengths` must be the same;

- The value `sum(baseline_lengths)` must be the same as the length of `pix_idx`.
"""
function baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm;
                          unseen=NaN)
    
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
    
    if(!ismissing(comm))
        noise_map = MPI.allreduce(partial_map, MPI.SUM, comm)
        hits = MPI.allreduce(partial_hits, MPI.SUM, comm)
    else 
        
        noise_map .= partial_map
        hits .= partial_hits
    end

    for i in eachindex(noise_map)
        if (hits[i] > 0)
            noise_map[i] /= hits[i]
        else
            noise_map[i] = unseen
        end
    end 
    
    noise_map
end


function applyz_and_sum(pix_idx, tod, baseline_lengths, num_of_pixels, comm;
                        unseen=NaN)
    @assert length(tod) == length(pix_idx)
        
    baselines_sum = zeros(eltype(tod), length(baseline_lengths))

    binned_map = tod2map_mpi(pix_idx, tod, num_of_pixels, comm, unseen=unseen)

    startidx = 1
    for i in eachindex(baseline_lengths)
        endidx = baseline_lengths[i] + startidx - 1

        for j in startidx:endidx
            baselines_sum[i] += tod[j] - binned_map[pix_idx[j]]
        end

        startidx += baseline_lengths[i]
    end
    
    baselines_sum
end


function applya(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm;
                unseen=NaN)
    @assert length(tod) == length(pix_idx)
    @assert length(baselines) == length(baseline_lengths)

    baselines_sum = zeros(eltype(baselines), length(baseline_lengths))
    total_sum = zero(eltype(baselines))

    binned_map = baseline2map_mpi(
        pix_idx,
        baselines,
        baseline_lengths,
        num_of_pixels,
        comm;
        unseen=unseen,
    )
    
    startidx = 1

    for i in eachindex(baseline_lengths)
        endidx = baseline_lengths[i] + startidx - 1
        
        for j in startidx:endidx
            baselines_sum[i] += baselines[i] - binned_map[pix_idx[j]]
        end

        startidx += baseline_lengths[i]
    end

    #needed to assure that sum(baselines)==0

    if(!ismissing(comm))   
        total_sum = MPI.allreduce([sum(baselines)], MPI.SUM, comm)[1]
    else
        total_sum = sum(baselines)
    end

    baselines_sum .+= total_sum 
end


function mpi_dot_prod(x, y, comm)

    @assert eltype(x) == eltype(y)
    result = zero(eltype(x))
    local_sum::eltype(x) = dot(x, y)

    if(!ismissing(comm))
        result = MPI.allreduce([local_sum], MPI.SUM, comm)[1]
    else
        result = local_sum
    end

    return result
end


function conj_grad(baselines_sum, pix_idx, tod, baseline_lengths,
                   num_of_pixels, rank, comm;
                   threshold=1e-9, max_iter=10000)
    
    T = eltype(tod)
    N = eltype(pix_idx)
    
    baselines = Array{T}(undef, length(baseline_lengths))
    r = Array{T}(undef, length(baseline_lengths))
    r_next = Array{T}(undef, length(baseline_lengths))
    p = Array{T}(undef, length(baseline_lengths))
    Ap = Array{T}(undef, length(baseline_lengths))
    
     
    baselines .= 0    #starting baselines 
    k = zero(N)       #number of iterations
    convergence_parameter = zero(T)
    rdotr = zero(T)
    rdotr_next = zero(T)

    best_convergence_parameter = zero(T)
    best_baselines = zeros(T, length(baseline_lengths))
    best_k = zero(N)
    
    r = baselines_sum - applya(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm)  #residual
    p .= r  
    rdotr = mpi_dot_prod(r, r, comm)
    best_convergence_parameter = sqrt(rdotr)

    if(rdotr == 0)
        return best_baselines
    end


    while true        

        Ap =  applya(p, pix_idx, tod, baseline_lengths, num_of_pixels, comm)

        rdotr = mpi_dot_prod(r, r, comm)
        pdotAp = mpi_dot_prod(p, Ap, comm)

        alpha = rdotr / pdotAp
        @. baselines += alpha * p  
        @. r_next = r - alpha * Ap
        
        rdotr_next = mpi_dot_prod(r_next, r_next, comm)
        
        convergence_parameter = sqrt(rdotr_next)
                
        if (convergence_parameter < best_convergence_parameter)
            best_convergence_parameter = convergence_parameter
            best_baselines .= baselines
            best_k = k
        end

        if convergence_parameter < threshold  break end
        if k  > max_iter   break end
        
        beta = rdotr_next/rdotr
        
        @. p = r_next + beta * p
        r .= r_next
        k += 1
    end

    return best_baselines
end


function destriped_map(baselines, pix_idx, tod,
                       baseline_lengths, num_of_pixels, comm;
                       unseen=NaN)
    
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, comm) - baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm)
    
end


@doc raw"""
    destripe(pix_idx, tod, num_of_pixels, baseline_lengths, rank, comm; threshold = 1e-9, max_iter = 10000) -> (pixels, baselines)

This MPI based function creates a map from a TOD and removes both 1/f
and white noise, using the destriping technique.

The parameters passed to the function have the following meaning:

- `pix_idx`: array containing the indices of the pixels visited by the
  instrument

- `tod`: the values measured by the polarimeters for each pixel
  (either I, Q, or U)

- `num_of_pixels`: the number of pixels in the map to be
  produced. This is used as an upper limit for the values in `pix_idx`

- `baseline_lengths`: number of samples in each baseline. The value
  `sum(baseline_lengths)` must be equal to `length(tod)`.

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

The function returns a 2-tuple containing the destriped map itself
(`Array{T,1}`, where `T` is the base type used for the TOD in the
argument`tod`) and an array containing the baselines.

Since it is not granted that the sequence of convergence parameters of
the conjugate gradient is monotonically decreasing, the code keeps the
lowest value of them and the corresponding array of baselines.  If the
loop ends because the maximum number of iterations has been reached,
this is the configuration that will be returned to the caller.

# Remarks

- The length of the arrays `pix_idx` and `tod` must be the same;

- `sum(baseline_lengths)` must be equal to `length(tod)`;

- If you are not using MPI, pass `missing` to the `comm` parameter.
"""
function destripe(pix_idx, tod, num_of_pixels, baseline_lengths, rank, comm;
                  threshold = 1e-9, max_iter = 10000, unseen=NaN)
    @assert sum(baseline_lengths) == length(tod)

    baselines_sum = applyz_and_sum(
        pix_idx,
        tod,
        baseline_lengths,
        num_of_pixels,
        comm,
        unseen=unseen,
    )
    
    baselines = conj_grad(
        baselines_sum,
        pix_idx,
        tod,
        baseline_lengths,
        num_of_pixels,
        rank,
        comm;
        threshold = threshold,
        max_iter = max_iter,
    )

    # once we have an estimate of the baselines, we can build the destriped map
    destr_map = destriped_map(
        baselines,
        pix_idx,
        tod,
        baseline_lengths,
        num_of_pixels,
        comm,
        unseen=unseen,
    )

    (destr_map, baselines)
end
