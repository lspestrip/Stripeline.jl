export tod2map_mpi, baseline2map_mpi, destripe

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


function applya(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
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




function conj_grad(baselines_sum, pix_idx, tod, baseline_lengths, num_of_pixels, comm; threshold = 1e-9, max_iter=10000)
    
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
    
    r = baselines_sum - applya(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm)  #residual

    p .= r  

    if(mpi_dot_prod(r, r, comm) == 0)
        return baselines
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

        if convergence_parameter < threshold  break end
        if k  > max_iter   break end
        
        beta = rdotr_next/rdotr
        
        @. p = r_next + beta * p
        r .= r_next
        k += 1
    end

    println("iteration number $k, residual = $convergence_parameter")

    return baselines
end



function destriped_map(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, comm) - baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels, comm)
end


"""
    destripe(pix_idx, tod, num_of_pixels, baseline_lengths, comm; threshold = 1e-9, max_iter = 10000) -> (pixels, baselines)

This MPI based function creates a map from a TOD and removes both 1/f and white noise, using the destriping technique. 

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map
-the array containg the length of each 1/f baseline
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

N.B.
* pix_idx and tod must be array of the same length and sum(baseline_lengths) must be equal to the length of `tod`.
* If you are not using MPI remember to initialize `comm` to `missing`.
"""
function destripe(pix_idx, tod, num_of_pixels, baseline_lengths, comm; threshold = 1e-9, max_iter = 10000, unseen=NaN)
    @assert sum(baseline_lengths) == length(tod)

    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_lengths, num_of_pixels, comm, unseen=unseen)
    baselines = conj_grad(baselines_sum, pix_idx, tod, baseline_lengths, num_of_pixels, comm; threshold = threshold, max_iter = max_iter)

    # once we have an estimate of the baselines, we can build the destriped map
    destr_map = destriped_map(baselines, pix_idx, tod, baseline_lengths, num_of_pixels, comm, unseen=unseen)

    #check that sum(baselines) = 0
    if(!ismissing(comm))   
        total_sum = MPI.allreduce([sum(baselines)], MPI.SUM, comm)[1]
    else
        total_sum = sum(baselines)
    end
    println("The sum of baselines is: $total_sum")

    (destr_map, baselines)
end