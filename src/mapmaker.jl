export tod2map, baseline2tod, applyz, applyz_and_sum, applya, destriped_map, destripe

using LinearMaps

import Healpix
import IterativeSolvers

#from TOD to binned map: it computes the mean value of the samples taken in each pixel
function tod2map(pix_idx, tod, num_of_pixels; unseen = NaN)
    binned_map = Array{Float64}(num_of_pixels)
    hits = zeros(Int, num_of_pixels)
    for i in eachindex(pix_idx)
        binned_map[pix_idx[i]] = binned_map[pix_idx[i]] + tod[i]
        hits[pix_idx[i]] = hits[pix_idx[i]] + 1
    end

    for i in eachindex(binned_map)
        if (hits[i] != 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end 
    
    binned_map
end

# Z = 1 - P*inv(M)*P'.  subtract the sky signal from the TOD
function applyz(pix_idx, tod, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    
    binned_map = tod2map(pix_idx, tod, num_of_pixels, unseen=unseen)
    tod - binned_map[pix_idx]
end

function applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    num_of_baselines = length(baseline_dim)
    baselines_sum = zeros(Float64, num_of_baselines)
    noise_tod = applyz(pix_idx, tod, num_of_pixels, unseen=unseen)

    count = 0
    for i in eachindex(baseline_dim)
        baselines_sum[i] = sum(noise_tod[(count + 1):(baseline_dim[i] + count)])
        count += baseline_dim[i]
    end
    
    baselines_sum
end

function baseline2tod(baselines, tod, baseline_dim)
    result = Array{Float64}(length(tod))

    count = 1
    for i in eachindex(baseline_dim)
        result[count:count + baseline_dim[i] - 1] = baselines[i]
        count += baseline_dim[i]
    end 

    result
end

function applya(a, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    baselines_tod = baseline2tod(a, tod, baseline_dim)
    applyz_and_sum(pix_idx, baselines_tod, baseline_dim, num_of_pixels, unseen=unseen)
end

function destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    baselines_tod = baseline2tod(baselines, tod, baseline_dim)
    tod2map(pix_idx, tod - baselines_tod, num_of_pixels, unseen=unseen)
end

function destripe(tod, pix_idx, baseline_dim, num_of_pixels; unseen=NaN)
    @assert sum(baseline_dim) == length(tod)

    A = LinearMap(x -> applya(x, pix_idx, tod, baseline_dim, num_of_pixels),
                  length(baseline_dim), length(baseline_dim),
                  issymmetric = true, ishermitian = true, isposdef = true)

    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)
    baselines = IterativeSolvers.cg(A, baselines_sum)

    # once we have an estimate of the baselines, we can build the destriped map
    pixels = destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)

    (pixels, baselines)
end
