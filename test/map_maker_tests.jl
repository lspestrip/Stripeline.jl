using Test

comm = missing

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
(pixels, baselines) = destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)
@test pixels ≈ true_map
@test baselines ≈ zeros(length(baselines))

# Signal only, but with an unseen pixel (the last one)
@test tod2map_mpi(pix_idx, tod, num_of_pixels + 1, comm, unseen = -1) ≈ [true_map; -1]

# Signal + noise
# This map only contains 1/f noise
@test tod2map_mpi(pix_idx, tod + oofnoise, num_of_pixels, comm) ≈ [4.25, 20.0, 10.0]

# Full process (destriping + map making)
(pixels, baselines) = destripe(pix_idx, tod + oofnoise, num_of_pixels, baseline_len, comm)
@test pixels ≈ true_map
@test baselines ≈ true_baselines
