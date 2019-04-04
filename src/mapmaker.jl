export data_properties_struct, get_data_properties, tod2map_mpi, baseline2map_mpi, destripe, baselines_covmat

using LinearAlgebra

try
    import MPI
catch
end 

@doc raw"""
This structure holds a number of parameters relative to the noise level and
the 1/f baselines measured or simulated for a certain polarimeter.

Field               | Type           | Meaning
:-----------------  |:-------------- |:----------------------------------------------------------------------------------
`pol_number`        | Int            | ID number of the polarimeter
`σ`                 | Float          | σ of white noise 
`num_of_samples`    | Int            | total number of samples
`num_of_baselines`  | Int            | total number of 1/f baselines in the measured or simulated TOD for this polarimeter
`baselines_lengths` | Array{Int}     | Array containing the length of each 1/f baseline for this polarimeter

# Example
data_properties(1, 0.002463, 500, 5, [100, 100, 100, 100, 100, 100]) 

says that we have 5 baselines simulated for polarimeter number 1 (in this case polarimeter STRIP31, with σ = 0.002463), each of length 100 samples.
"""
struct data_properties_struct
    polarimeter::Int
    σ::Float64
    number_of_samples::Int
    number_of_baselines::Int
    baselines_lengths::Array{Int}
end


@doc raw"""
    function get_data_properties(detector_number, σ_k, num_of_baselines, num_of_samples) -> data_properties
    
This function can be used to build the object `data properties` needed for desriping and calculation of baselines_covmat.

It requires in input 4 arrays containing:
- the ID number of the polarimeters that the current rank will simulate
- the white noise σ for each polarimeter that the current rank will simulate
- the number of 1/f baselines for each polarimeter that the current rank will simulate
- the total number of samples for each polarimeter that the current rank will simulate

N.B. the ID numbers, the number of 1/f baselines and the total number of samples can be obtained using the function `get_chunk_properties`.

It returns an array of `data_properties_struct`, of length equal to the number of polarimeters simulated by current rank.

"""
function get_data_properties(detector_number, σ_k, num_of_baselines, num_of_samples)
    data_properties  = Array{data_properties_struct}(undef, length(detector_number))
    for i in 1:length(detector_number)
        baselines_lengths  = repeat([Int64(num_of_samples[i]/num_of_baselines[i])], num_of_baselines[i])  #we suppose baselines of equal lengths
        data_properties[i]  = data_properties_struct(detector_number[i], σ_k[detector_number[i]], num_of_samples[i], num_of_baselines[i], baselines_lengths)
    end
    data_properties
end



@doc raw"""
    tod2map(pix_idx, tod, num_of_pixels, comm) -> binned_map
    tod2map(pix_idx, tod, num_of_pixels, data_properties, comm) -> binned_map

This function creates a binned map from the time-ordered data kept in
the array `tod`, assuming that each sample is observing the pixel
whose index is in `pix_idx`. The parameter `num_of_pixels` contains
the number of pixels in the Healpix map, and it is used as an upper
bound for the values in `pix_idx`. The parameter `comm` must be a MPI
communicator, or `missing` if you are not using MPI.
This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with MPI.allreduce.
The function returns an array containing the binned map.

If Array of structures `data_properties_struct`is passed to the function, the output binned map will be a weighted binned map.
Each sample will be weighted according to the inverse white noise variance σ^2 of the corrisponding polarimeter.
In this way, the less noisy polarimeters will count more in the estimation of the map.

# Requirements
- The length of the arrays `pix_idx` and `tod` must be the same

"""
tod2map_mpi


function tod2map_mpi(pix_idx, tod, num_of_pixels, comm; unseen = NaN)   
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

function tod2map_mpi(pix_idx, tod, num_of_pixels, data_properties, comm; unseen = NaN)   
    T = eltype(tod)
    
    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(T, num_of_pixels)
    binned_map = zeros(T, num_of_pixels)
    hits = zeros(T, num_of_pixels )

    start_idx = 1

    for j in eachindex(data_properties)                 #loop on detectors
        end_idx = start_idx + data_properties[j].number_of_samples-1
        for i in start_idx:end_idx             #loop on samples 
            partial_map[pix_idx[i]] += tod[i] * 1/(data_properties[j].σ)^2
            partial_hits[pix_idx[i]] += 1/(data_properties[j].σ)^2
        end
        start_idx += data_properties[j].number_of_samples
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
    baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm)-> noise_map
    baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties, comm) -> noise_map

This is a MPI based function: each MPI process computes a map from its
available data.  All partial maps are then combined together with
MPI.allreduce.
The function returns an array containing the binned map.

If Array of structures `data_properties_struct`is passed to the function (instead of `baseline_lengths`) the output binned map will be a weighted binned map.
Each sample will be weighted according to the inverse white noise variance σ^2 of the corrisponding polarimeter.
In this way, the less noisy polarimeters will count more in the estimation of the map.

# Requirements
- The length of `baselines` and `baseline_lengths` must be the same;
- The value `sum(baseline_lengths)` must be the same as the length of `pix_idx`.
"""
baseline2map_mpi


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

function baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties, num_of_baselines, comm; unseen=NaN)
    
    T = eltype(baselines)
    
    partial_map = zeros(T, num_of_pixels)
    partial_hits = zeros(T, num_of_pixels)
    noise_map = zeros(T, num_of_pixels)
    hits = zeros(T, num_of_pixels)

    startidx = 1
    baseline_idx = 1 

    for l in eachindex(data_properties)  #loop on detectors
        for i in eachindex(data_properties[l].baselines_lengths)
            endidx = data_properties[l].baselines_lengths[i] + startidx - 1        
        
            for j in startidx:endidx
                partial_map[pix_idx[j]] += baselines[baseline_idx] * 1/(data_properties[l].σ)^2
                partial_hits[pix_idx[j]] += 1/(data_properties[l].σ)^2
            end
            startidx += data_properties[l].baselines_lengths[i]
            baseline_idx += 1
        end 
    end



function applyz_and_sum(pix_idx, tod, num_of_pixels, data_properties, num_of_baselines, comm; unseen=NaN)

    @assert length(tod) == length(pix_idx)
        
    baselines_sum = zeros(eltype(tod), num_of_baselines)

    binned_map = tod2map_mpi(
        pix_idx, 
        tod, 
        num_of_pixels, 
        data_properties, 
        comm, 
        unseen=unseen)

    baseline_idx = 1
    startidx = 1
    for l in eachindex(data_properties)  #loop on detectors
      
        for i in eachindex(data_properties[l].baselines_lengths)
            endidx = data_properties[l].baselines_lengths[i] + startidx - 1

            for j in startidx:endidx
                baselines_sum[baseline_idx] += (tod[j] - binned_map[pix_idx[j]]) * 1/(data_properties[l].σ)^2
            end

            startidx += data_properties[l].baselines_lengths[i]
            baseline_idx += 1
        end
    end
    
    baselines_sum
end


function applya(baselines, pix_idx, num_of_baselines, num_of_pixels, data_properties, comm; unseen=NaN)
    @assert length(baselines) == num_of_baselines
    
    baselines_sum = zeros(eltype(baselines), num_of_baselines)
    total_sum = zero(eltype(baselines))

    binned_map = baseline2map_mpi(
        pix_idx, 
        baselines, 
        num_of_pixels, 
        data_properties, 
        num_of_baselines, 
        comm;
        unseen=unseen)

    startidx = 1
    baseline_idx = 1

    for l in eachindex(data_properties)
      
        for i in eachindex(data_properties[l].baselines_lengths)
            endidx = data_properties[l].baselines_lengths[i] + startidx - 1
        
            for j in startidx:endidx
                baselines_sum[baseline_idx] += (baselines[baseline_idx] - binned_map[pix_idx[j]]) * 1/(data_properties[l].σ)^2
            end

            startidx += data_properties[l].baselines_lengths[i]
            baseline_idx += 1
        end
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



function conj_grad(baselines_sum, pix_idx, tod,
                   num_of_pixels, data_properties, 
                   num_of_baselines, rank, comm; 
                   threshold = 1e-9, max_iter=10000)

    T = eltype(tod)
    N = eltype(pix_idx)
    
    baselines = Array{T}(undef, num_of_baselines)
    r = Array{T}(undef, num_of_baselines)
    r_next = Array{T}(undef, num_of_baselines)
    p = Array{T}(undef, num_of_baselines)
    Ap = Array{T}(undef,num_of_baselines)
    
     
    baselines .= 0    #starting baselines 
    k = zero(N)       #number of iterations
    convergence_parameter = zero(T)
    rdotr = zero(T)
    rdotr_next = zero(T)

    best_convergence_parameter = zero(T)
    best_baselines = zeros(T, num_of_baselines)
    best_k = zero(N)

    r = baselines_sum - applya(baselines, pix_idx, num_of_baselines, num_of_pixels, data_properties, comm)  #residual
    p .= r  
    rdotr = mpi_dot_prod(r, r, comm)
    best_convergence_parameter = sqrt(rdotr)

    if(rdotr == 0)
        return best_baselines
    end


    while true        

        Ap =  applya(p, pix_idx,  num_of_baselines, num_of_pixels, data_properties, comm)

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


function destriped_map(baselines, pix_idx, tod, data_properties, num_of_pixels, num_of_baselines, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, data_properties, comm) - baseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties, num_of_baselines, comm)
end




@doc raw"""
    destripe(pix_idx, tod, num_of_pixels, data_properties, rank, comm; threshold = 1e-9, max_iter = 10000) -> (pixels, baselines)

This MPI based function creates a map from a TOD and removes both 1/f
and white noise, using the destriping technique.

The parameters passed to the function have the following meaning:

- `pix_idx`: array containing the indices of the pixels visited by the instrument

- `tod`: the values measured by the polarimeters for each pixel
  (either I, Q, or U)

- `num_of_pixels`: the number of pixels in the map to be
   produced. This is used as an upper limit for the values in `pix_idx`

- `data_properties`: an array of structures `data_properties_struct`, 
   holding information on each simulated polarimeter noise level, 
   number and length of 1/f baselines and total number of samples.
   It can be obtained by using function `get_data_properties`.

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
- If you are not using MPI, pass `missing` to the `comm` parameter.
"""
function destripe(pix_idx, tod, num_of_pixels, data_properties, rank, comm;
                 threshold = 1e-9, max_iter = 10000, unseen=NaN)

    num_of_baselines = 0
    for i in 1:length(data_properties) num_of_baselines  += data_properties[i].number_of_baselines end

    baselines_sum = applyz_and_sum(
            pix_idx, 
            tod,
            num_of_pixels, 
            data_properties, 
            num_of_baselines, 
            comm, 
            unseen=unseen)

    baselines = conj_grad(
            baselines_sum, 
            pix_idx, 
            tod, 
            num_of_pixels, 
            data_properties, 
            num_of_baselines, 
            rank, 
            comm; 
            threshold = threshold, 
            max_iter = max_iter)

    # once we have an estimate of the baselines, we can build the destriped map
    destr_map = destriped_map(
            baselines,
            pix_idx, 
            tod, 
            data_properties, 
            num_of_pixels,
            num_of_baselines, 
            comm, 
            unseen=unseen)

    #check that sum(baselines) = 0
    if(!ismissing(comm))   
        total_sum = MPI.allreduce([sum(baselines)], MPI.SUM, comm)[1]
    else
        total_sum = sum(baselines)
    end

    if rank==0
        println("The sum of baselines is: $total_sum")
    end

    (destr_map, baselines)
end


@doc raw"""
    baselines_covmat(polarimeters, σ_k, baseline_length_s, fsamp_hz, total_time) -> covariance_matrix

This function produces the covariance matrix of the 1/f baselines computed by the destriper.
As an approximation, we ignore nondiagonal terms (i.e. the correlations between different baselines).
The function output is therefore an array corresponding to the covariance matrix diagonal.
Another assumption is that the white noise variance stays constant over a given baseline.

The baseline error is computed according to equation 29 in https://arxiv.org/abs/0904.3623

The parameters passed to the function have the following meaning:
    -`polarimeters`: the array of polarimeters ID numbers 
    -`σ_k`: the array of corresponding white noise σ (in K)
    -`baseline_length_s` = the length (in s) of each 1/f baseline
    -`fsamp_hz` = the sampling frequency (in Hz)
    -`total_time` = the duration (in s) of the observation

The function firstly computes an Array of structures `data_properties_struct`,
which matches the white noise variance σ with the number of 1/f baselines and their lengths, for each polarimeter.
Then it computes the covariance matrix according to these matches.     
"""       
function baselines_covmat(polarimeters, σ_k, baseline_length_s, fsamp_hz, total_time)

    ALL_data_properties = Array{data_properties_struct}(undef, length(polarimeters))
    baselines_per_pol =  Int64(total_time/baseline_length_s)
    for i in 1:length(polarimeters)
        baselines_lengths  = repeat([baseline_length_s*fsamp_hz], baselines_per_pol)
        ALL_data_properties[i]  = data_properties_struct(polarimeters[i], σ_k[i], sum(baselines_lengths), baselines_per_pol, baselines_lengths)
    end

    covariance_matrix = []
    for i in 1:length(ALL_data_properties) 
        partial_covmat = Array{Float64}(undef,length(ALL_data_properties[i].baselines_lengths))
        
        for j in 1:length(ALL_data_properties[i].baselines_lengths) 
            σ = ALL_data_properties[i].σ
            this_baseline_length = ALL_data_properties[i].baselines_lengths[j]  
            
            partial_covmat[j] = σ^2/this_baseline_length    
        end
        covariance_matrix = append!(covariance_matrix, partial_covmat)
    end
    return covariance_matrix
end


   