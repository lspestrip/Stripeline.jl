using Test

comm = missing
rank = 0

# Index of the pixels: this fixes the length of the TOD as well
pix_idx = [1, 1, 2, 1, 2, 3, 3, 1, 3, 2]
# The sky map
true_map = [4., 20., 10.]
num_of_pixels = length(true_map)

# This is the realization of the 1/f noise we're going to use later
true_baselines = [1., -2, 1.]
baseline_len = [3, 3, 4]

# This TOD does not include 1/f noise
tod = true_map[pix_idx]

# This TOD includes 1/f noise only. This must match the definition of "true_baselines" and "baseline_len" above
oofnoise = [1., 1., 1., -2., -2., -2., 1., 1., 1., 1.]


# Signal only
@test tod2map_mpi(pix_idx, tod, num_of_pixels, comm) ≈ true_map
(pixels, baselines) = destripe(pix_idx, tod, num_of_pixels, baseline_len, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ zeros(length(baselines))

# Signal only, but with an unseen pixel (the last one)
@test tod2map_mpi(pix_idx, tod, num_of_pixels + 1, comm, unseen = -1) ≈ [true_map; -1]

# Signal + noise
# This map only contains 1/f noise
@test tod2map_mpi(pix_idx, tod + oofnoise, num_of_pixels, comm) ≈ [4.25, 20.0, 10.0]

# Full process (destriping + map making)
(pixels, baselines) = destripe(pix_idx, tod + oofnoise, num_of_pixels, baseline_len, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ true_baselines


#check divergences
pix_idx = repeat([1, 1, 2, 1, 2, 3, 3, 1, 3, 2],1000)
baseline_len = repeat([3,3,4],1000)
ooftod = repeat([1., 1., 1., -2., -2., -2., 1., 1., 1., 1.],1000)
skytod = repeat([3.86],10000)
true_map = [3.86, 3.86, 3.86]
true_baselines = repeat([1,-2, 1],1000)

(pixels, baselines) = destripe(pix_idx, skytod+ooftod, num_of_pixels, baseline_len, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ true_baselines


#Test covariance matrix of baselines computation
polarimeters = [8, 48, 67] 
σ_k = [0.0020967137443360585, 0.003495923350914105, 0.0016551546163219948]
baseline_length_s = 10
fsamp_hz = 10
total_time = 20
covmat = baselines_covmat(polarimeters, σ_k, baseline_length_s, fsamp_hz, total_time)

@test covmat ≈ [4.396208525687735e-8 , 4.396208525687735e-8 , 1.2221480075466504e-7, 1.2221480075466504e-7, 2.73953680393201e-8, 2.73953680393201e-8]