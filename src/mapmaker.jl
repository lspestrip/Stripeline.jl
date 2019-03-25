export tod2map_mpi, baseline2map_mpi, destripe, baselines_covmat

using LinearAlgebra

try
    import MPI
catch
end 

"""
    tod2map(pix_idx, tod, num_of_pixels) -> binned_map
This function creates a binned map from a TOD, removing white noise.
This is a MPI based function: each MPI process computes a map from its available data.
All partial maps are then combined together with MPI.allreduce.

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map
-the MPI communicator

N.B.
* pix_idx and tod must be array of the same length.
* If you are not using MPI remember to initialize `comm` to `missing`

"""
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


function baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    
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


function applyz_and_sum(pix_idx, tod, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
        
    baselines_sum = zeros(eltype(tod), length(baseline_lengths))

    binned_map = tod2map_mpi(pix_idx, tod, num_of_pixels, comm, unseen=unseen)

    startidx = 1
    for i in eachindex(baseline_lengths)
        endidx = baseline_lengths[i] + startidx - 1

        # The inner for is equivalent to
        #
        #   baselines_sum[i] += sum(tod[startidx:endidx] - binned_map[pix_idx[startidx:endidx]])
        #
        # but roundoff errors are reduced
        for j in startidx:endidx
            baselines_sum[i] += tod[j] - binned_map[pix_idx[j]]
        end

        startidx += baseline_lengths[i]
    end
    
    baselines_sum
end


function applya(baselines, pix_idx, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    @assert length(baselines) == length(baseline_lengths)

    baselines_sum = zeros(eltype(baselines), length(baseline_lengths))
    total_sum = zero(eltype(baselines))

    binned_map = baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    
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




function conj_grad(baselines_sum, pix_idx, tod, baseline_lengths, num_of_pixels, rank, comm; threshold = 1e-9, max_iter=10000)
    
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
    
    r = baselines_sum - applya(baselines, pix_idx, baseline_lengths, num_of_pixels, comm)  #residual
    p .= r  
    rdotr = mpi_dot_prod(r, r, comm)
    best_convergence_parameter = sqrt(rdotr)

    if(rdotr == 0)
        return best_baselines
    end


    while true        

        Ap =  applya(p, pix_idx, baseline_lengths, num_of_pixels, comm)

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

    if rank==0
        println("Last iteration number $k, Last residual = $convergence_parameter")
        println("BEST iteration number $best_k, BEST residual = $best_convergence_parameter")
    end

    return best_baselines
end



function destriped_map(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, comm) - baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm)
end


"""
    destripe(pix_idx, tod, num_of_pixels, baseline_lengths, rank, comm; threshold = 1e-9, max_iter = 10000) -> (pixels, baselines)

This MPI based function creates a map from a TOD and removes both 1/f and white noise, using the destriping technique. 

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map
-the array containg the length of each 1/f baseline
-the MPI rank number
-the MPI communicator

and, as optional arguments:
-the conjugate gradient threshold: when the residual error of the iteration goes
below this value, the iteration stops. The smaller the value,
the more accurate the solution. 
Default = 1e-09
-the maximum number of iterations: if the CG algorithm does not converge after
this number of steps, quit the iteration.
Default = 10000

It returns a tuple containing the destriped map itself (Array{Float64,1}) and the estimated array of 1/f baselines.

Since it is not granted that the sequence of convergence parameters of the conjugate gradient is
monotonically decreasing, the code keeps the lowest value of them and the corresponding array of baselines.
If the loop ends because the maximum number of iterations has been reached, this is the configuration that will be returned to
the caller.

N.B.
* pix_idx and tod must be array of the same length and sum(baseline_lengths) must be equal to the length of `tod`.
* If you are not using MPI remember to initialize `comm` to `missing`.
"""
function destripe(pix_idx, tod, num_of_pixels, baseline_lengths, rank, comm; threshold = 1e-9, max_iter = 10000, unseen=NaN)
    @assert sum(baseline_lengths) == length(tod)

    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_lengths, num_of_pixels, comm, unseen=unseen)
    baselines = conj_grad(baselines_sum, pix_idx, tod, baseline_lengths, num_of_pixels, rank, comm; threshold = threshold, max_iter = max_iter)

    # once we have an estimate of the baselines, we can build the destriped map
    destr_map = destriped_map(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm, unseen=unseen)

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




"""
This structure holds a number of parameters relative to the 1/f baselines measured or simulated for a certain polarimeter.
It is useful for the computation of the covariance matrix of baselines.

Field               | Type           | Meaning
:-----------------  |:-------------- |:----------------------------------------------------------------------------------
`pol_number`        | Int            | ID number of the polarimeter
`σ`                 | Float          | σ of white noise 
`num_of_baselines`  | Int            | total number of 1/f baselines in the measured or simulated TOD for this polarimeter
`baselines_lengths` | Array{Int}     | Array containing the length of each 1/f baseline for this polarimeter

# Example
match_σ_baselines(31, 0.002463, 5, [100, 100, 100, 100, 100, 100]) 

says that we have 5 baselines simulated for polarimters number 31 (with σ = 0.002463), each of length 100 samples.
"""
struct match_σ_baselines
    polarimeter::Int
    σ::Float64
    number_of_baselines::Int
    baselines_lengths::Array{Int}
end 



"""
baselines_covmat(polarimeters, σ_k, baseline_length_s, fsamp_hz, total_time) -> covariance_matrix

This function produces the covariance matrix of the 1/f baselines computed by the destriper.
As an approximation, we ignore nondiagonal terms (i.e. the correlations between different baselines).
The function output is therefore an array corresponding to the covariance matrix diagonal.
Another assumption is that the white noise variance stays constant over a given baseline.

The baseline error is computed according to equation 29 in https://arxiv.org/abs/0904.3623

The function requires in input:
-the array of polarimeters ID numbers 
-the array of corresponding white noise σ (in K)
-the length (in s) of each 1/f baseline
-the sampling frequency (in Hz)
-the duration (in s) of the observation

The function firstly computes an Array of structures `match_σ_baselines`, which matches the white noise variance σ with the number of 1/f baselines
and their lengths, for each polarimeter.
Then it computes the covariance matrix according to these matches.     
"""       
function baselines_covmat(polarimeters, σ_k, baseline_length_s, fsamp_hz, total_time)

    matches_σ_baselines  = Array{match_σ_baselines}(undef, length(polarimeters))
    baselines_per_pol =  Int64(total_time/baseline_length_s)
    for i in 1:length(polarimeters)
        baselines_lengths  = repeat([baseline_length_s*fsamp_hz], baselines_per_pol)
        matches_σ_baselines[i]  = match_σ_baselines(polarimeters[i], σ_k[i], baselines_per_pol, baselines_lengths)
    end
        
    covariance_matrix = []
    for i in 1:length(matches_σ_baselines) 
        partial_covmat = Array{Float64}(undef,length(matches_σ_baselines[i].baselines_lengths))
        
        for j in 1:length(matches_σ_baselines[i].baselines_lengths) 
            σ = matches_σ_baselines[i].σ
            this_baseline_length = matches_σ_baselines[i].baselines_lengths[j]  
            
            partial_covmat[j] = σ^2/this_baseline_length    
        end
        covariance_matrix = append!(covariance_matrix, partial_covmat)
    end
    return covariance_matrix
end
