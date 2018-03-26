using Base.Test

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

# This TOD only contains 1/f noise
oofnoise = baseline2tod(true_baselines, tod, baseline_len)

# This must match the definition of "true_baselines" and "baseline_len" above
@test oofnoise == [1., 1., 1., -2., -2., -2., 1., 1., 1., 1.]

# Signal only
@test tod2map(pix_idx, tod, num_of_pixels) ≈ true_map
@test applyz_and_sum(pix_idx, tod, baseline_len, num_of_pixels) == zeros(length(true_baselines))
@test applya(true_baselines, pix_idx, tod, baseline_len, num_of_pixels) ≈ [2.5, -6.25, 3.75]
(pixels, baselines) = destripe(tod, pix_idx, baseline_len, num_of_pixels)
@test pixels ≈ true_map
@test baselines ≈ zeros(length(baselines))

# Signal only, but with an unseen pixel (the last one)
@test tod2map(pix_idx, tod, num_of_pixels + 1, unseen = -1) ≈ [true_map; -1]
@test applyz_and_sum(pix_idx, tod, baseline_len, num_of_pixels + 1, unseen = -1) == zeros(length(true_baselines))
@test applya(true_baselines, pix_idx, tod, baseline_len, num_of_pixels + 1, unseen = -1) ≈ [2.5, -6.25, 3.75]

# Signal + noise
@test tod2map(pix_idx, tod + oofnoise, num_of_pixels) ≈ [4.25, 20.0, 10.0]
@test applyz_and_sum(pix_idx, tod + oofnoise, baseline_len, num_of_pixels) ≈ [2.5, -6.25, 3.75]
@test applya(true_baselines, pix_idx, tod + oofnoise, baseline_len, num_of_pixels) ≈ [2.5, -6.25, 3.75]
@test destriped_map(true_baselines, pix_idx, tod + oofnoise, baseline_len, num_of_pixels) ≈ true_map

# Full process (destriping + map making)
(pixels, baselines) = destripe(tod + oofnoise, pix_idx, baseline_len, num_of_pixels)
@test pixels ≈ true_map
@test baselines ≈ true_baselines