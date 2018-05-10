export tod2map_mpi, destripe


"""
    tod2map(pix_idx, tod, num_of_pixels) -> binned_map
This function create a binned map from a TOD, removing white noise.
This is a MPI based function: each MPI process computes a map from its available data.
All partial maps are then combined together with MPI.allreduce.

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map

pix_idx and tod must be array of the same length.

"""
function tod2map_mpi(pix_idx, tod, num_of_pixels; unseen = NaN)   

    partial_map = zeros(Float64,num_of_pixels)
    partial_hits = zeros(Int, num_of_pixels)

    for i in eachindex(pix_idx)
        partial_map[pix_idx[i]] = partial_map[pix_idx[i]] + tod[i]
        partial_hits[pix_idx[i]] = partial_hits[pix_idx[i]] + 1
    end

    binned_map = MPI.allreduce(partial_map, MPI.SUM, comm)
    hits = MPI.allreduce(partial_hits, MPI.SUM, comm)

    for i in eachindex(binned_map)
        if (hits[i] != 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end 
    
    binned_map
end


function baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels; unseen=NaN)

    partial_map = zeros(Float64, num_of_pixels)
    partial_hits = zeros(Int, num_of_pixels)
    startidx = 1
    
    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1        
        
        for j in startidx:endidx
            partial_map[pix_idx[j]] = partial_map[pix_idx[j]] +  baselines[i]
            partial_hits[pix_idx[j]] = partial_hits[pix_idx[j]] + 1
        end
        startidx += baseline_dim[i]
    end 

    noise_map = MPI.allreduce(partial_map, MPI.SUM, comm)
    hits = MPI.allreduce(partial_hits, MPI.SUM, comm)
    
    for i in eachindex(noise_map)
        if (hits[i] != 0)
            noise_map[i] = noise_map[i] / hits[i]
        else
            noise_map[i] = unseen
        end
    end 
    
    noise_map
end


function applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = tod2map_mpi(pix_idx, tod, num_of_pixels, unseen=unseen)

    startidx = 1
    baselines_sum = zeros(Float64, length(baseline_dim))

    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1

        # The inner for is equivalent to
        #
        #   baselines_sum[i] += sum(tod[startidx:endidx] - binned_map[pix_idx[startidx:endidx]])
        #
        # but roundoff errors are reduced
        for j in startidx:endidx
            baselines_sum[i] += tod[j] - binned_map[pix_idx[j]]
        end

        startidx += baseline_dim[i]
    end
    
    baselines_sum
end


function applya(baselines, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels; unseen=NaN)
    startidx = 1
    baselines_sum = zeros(Float64, length(baseline_dim))

    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1
        
        for j in startidx:endidx
            baselines_sum[i] += baselines[i] - binned_map[pix_idx[j]]
        end

        startidx += baseline_dim[i]
    end
    
    baselines_sum
end


function mpi_dot_prod(x, y)
    local_sum = dot(x, y)
    total = MPI.allreduce([local_sum], MPI.SUM, comm)[1]
    return total
end


function conj_grad(baselines_sum, start_baselines, pix_idx, tod, baseline_dim, num_of_pixels; threshold = 1e-9, max_iter=10000)
    #initialization 
    k = 0
    x = start_baselines
    r = baselines_sum - applya(start_baselines, pix_idx, tod, baseline_dim, num_of_pixels)  #residual
    p = r
    old_r_dot  = mpi_dot_prod(r, r) 

    while true
        k+=1  #iteration number
        if k >= max_iter 
            break
        end

        Ap =  applya(p, pix_idx, tod, baseline_dim, num_of_pixels)
        alpha = old_r_dot/mpi_dot_prod(p, Ap)      
        x = x + alpha*p
        r = r - alpha*Ap
        new_r_dot = mpi_dot_prod(r, r)

        if sqrt(new_r_dot) < threshold
            break
        end

        beta = new_r_dot/old_r_dot
        p = r + beta*p
        old_r_dot = new_r_dot
    end
    return x
end



function destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels) - baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels)
end


"""
    destripe(pix_idx, tod, num_of_pixels::Integer, baseline_dim) -> (pixels, baselines)

This MPI based function create a map from a TOD and removes both 1/f and white noise, using the destriping technique. 

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map
-the array containg the dimension of each 1/f baseline.

pix_idx and tod must be array of the same length and sum(baseline_dim) must be equal to the length of `tod`.

It returns a tuple containing the destriped map itself and the estimated array of 1/f baselines.
    
"""
function destripe(pix_idx, tod, num_of_pixels, baseline_dim; unseen=NaN)
    @assert sum(baseline_dim) == length(tod)

    start_baselines = zeros(Float64, length(baseline_dim))
    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)
    baselines = conj_grad(baselines_sum, start_baselines, pix_idx, tod, baseline_dim, num_of_pixels)

    # once we have an estimate of the baselines, we can build the destriped map
    pixels = destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)

    (pixels, baselines)
end