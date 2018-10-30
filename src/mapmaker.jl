export tod2map_mpi, destripe

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
    partial_map = zeros(Float64,num_of_pixels)
    partial_hits = zeros(Int, num_of_pixels)

    for i in eachindex(pix_idx)
        partial_map[pix_idx[i]] += tod[i]
        partial_hits[pix_idx[i]] += 1
    end

    if(!ismissing(comm))
        binned_map::Array{Float64, 1}= MPI.allreduce(partial_map, MPI.SUM, comm)
        hits::Array{Int, 1} = MPI.allreduce(partial_hits, MPI.SUM, comm)
    else
        binned_map = partial_map
        hits = partial_hits
    end

    @inbounds for i in eachindex(binned_map)
        if (hits[i] != 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end

    binned_map
end


function baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels, comm; unseen=NaN)

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
    
    if(!ismissing(comm))
        noise_map::Array{Float64, 1}= MPI.allreduce(partial_map, MPI.SUM, comm)
        hits::Array{Int, 1} = MPI.allreduce(partial_hits, MPI.SUM, comm)
    else 
        
        noise_map = partial_map
        hits = partial_hits
    end

    for i in eachindex(noise_map)
        if (hits[i] != 0)
            noise_map[i] = noise_map[i] / hits[i]
        else
            noise_map[i] = unseen
        end
    end 
    
    noise_map
end


function applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = tod2map_mpi(pix_idx, tod, num_of_pixels, comm, unseen=unseen)

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


function applya(baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels, comm; unseen=NaN)
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


function mpi_dot_prod(x, y, comm)
    local_sum = dot(x, y)

    if(!ismissing(comm))
        total::Float64 = MPI.allreduce([local_sum], MPI.SUM, comm)[1]
    else
        total = local_sum
    end

    return total
end


function conj_grad(baselines_sum, start_baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm; threshold = 1e-9, max_iter=10000) 
    k = 0
    x = start_baselines
    r = baselines_sum .- applya(start_baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm)  #residual
    p = r
    old_r_dot  = mpi_dot_prod(r, r, comm)
    
    if old_r_dot == 0
        return x
    end

    while true
        k+=1  #iteration number

        if k >= max_iter 
            break
        end

        Ap =  applya(p, pix_idx, tod, baseline_dim, num_of_pixels, comm)
        alpha = old_r_dot/mpi_dot_prod(p, Ap, comm) 
        x .+= alpha .* p
        r .-= alpha .* Ap
        new_r_dot = mpi_dot_prod(r, r, comm)

        if sqrt(new_r_dot) < threshold
            break
        end

        beta = new_r_dot/old_r_dot
        p = r .+ beta .* p
        old_r_dot = new_r_dot
    end
    return x
end



function destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map_mpi(pix_idx, tod, num_of_pixels, comm) - baseline2map_mpi(pix_idx, baselines, baseline_dim, num_of_pixels, comm)
end


"""
    destripe(pix_idx, tod, num_of_pixels::Integer, baseline_dim) -> (pixels, baselines)

This MPI based function creates a map from a TOD and removes both 1/f and white noise, using the destriping technique. 

It requires in input:
-the array of pointed pixels
-the TOD
-the desired number of pixels of the output map
-the array containg the dimension of each 1/f baseline
-the MPI communicator

It returns a tuple containing the destriped map itself (Array{Float64,1}) and the estimated array of 1/f baselines.

N.B.
* pix_idx and tod must be array of the same length and sum(baseline_dim) must be equal to the length of `tod`.
* If you are not using MPI remember to initialize `comm` to `missing`.
"""
function destripe(pix_idx, tod, num_of_pixels, baseline_dim, comm; unseen=NaN)
    @assert sum(baseline_dim) == length(tod)

    start_baselines = zeros(Float64, length(baseline_dim))
    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels, comm, unseen=unseen)
    baselines = conj_grad(baselines_sum, start_baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm)

    # once we have an estimate of the baselines, we can build the destriped map
    pixels = destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels, comm, unseen=unseen)

    (pixels, baselines)
end




