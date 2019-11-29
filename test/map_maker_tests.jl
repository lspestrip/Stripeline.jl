using Test
using Healpix
using Statistics 

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

data_properties = [TodNoiseProperties(pol = 1, rms = 1., baselines = baseline_len)]

# Signal only
@test tod2map_mpi(pix_idx, tod, num_of_pixels, data_properties, comm) ≈ true_map
(pixels, baselines) = destripe(pix_idx, tod, num_of_pixels, data_properties, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ zeros(length(baselines))

# Signal only, but with an unseen pixel (the last one)
@test tod2map_mpi(pix_idx, tod, num_of_pixels + 1, data_properties, comm, unseen = -1) ≈ [true_map; -1]

# Signal + noise
# This map only contains 1/f noise
@test tod2map_mpi(pix_idx, tod + oofnoise, num_of_pixels, data_properties, comm) ≈ [4.25, 20.0, 10.0]

# Full process (destriping + map making)
(pixels, baselines) = destripe(pix_idx, tod + oofnoise, num_of_pixels, data_properties, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ true_baselines


#check divergences
pix_idx = repeat([1, 1, 2, 1, 2, 3, 3, 1, 3, 2], 1000)
baseline_len = repeat([3,3,4], 1000)
ooftod = repeat([1., 1., 1., -2., -2., -2., 1., 1., 1., 1.], 1000)
skytod = repeat([3.86], 10000)
true_map = [3.86, 3.86, 3.86]
true_baselines = repeat([1,-2, 1], 1000)

data_properties = [TodNoiseProperties(pol = 1, rms = 1., baselines = baseline_len)]

(pixels, baselines) = destripe(pix_idx, skytod + ooftod, num_of_pixels, data_properties, rank, comm)
@test pixels ≈ true_map
@test baselines ≈ true_baselines

# Test covariance matrix of baselines computation

polarimeters = [8, 48, 67] 
sigma_k = [0.0020967137443360585, 0.003495923350914105, 0.0016551546163219948]
baseline_length_s = 10
fsamp_hz = 10
total_time = 20
covmat = baselines_covmat(polarimeters, sigma_k, baseline_length_s, fsamp_hz, total_time)

@test covmat ≈ [
    4.3962085256877350e-8,
    4.3962085256877350e-8,
    1.2221480075466504e-7,
    1.2221480075466504e-7,
    2.7395368039320100e-8,
    2.7395368039320100e-8,
]

# Test weighting destriper

comm = missing
commsize = 1
rank = 0
total_time = 60
fsamp_hz = 10
baseline_length_s = 10
tau_s = 1. / fsamp_hz
NSIDE = 32
resol = Healpix.Resolution(NSIDE) #desired resolution for output map
num_of_pixels = resol.numOfPixels

# 2 polarimeters 
num_of_polarimeters = 2
sigma_k =  [0.001, 1000]
fknee_hz = [0.01, 0.01]
slope = [1,-1]
baselines_per_process = Int64(total_time / baseline_length_s) * num_of_polarimeters

chunks = Any[Any[datachunk(1, 1, 6, 6), datachunk(2, 1, 6, 6)]]
data_properties = [
    TodNoiseProperties(pol = 1, rms = 0.001, baselines = [100, 100, 100, 100, 100, 100]), 
    TodNoiseProperties(pol = 2, rms = 1000.0, baselines = [100, 100, 100, 100, 100, 100]),
]

pix_idx = Int64[]

for i in 1:2  #loop on detectors
    times = 0.0:tau_s:total_time - tau_s
    (dirs, ψ) = genpointings([0.,0.,1.], times; latitude_deg = 28.29) do time_s
        return (0.0, deg2rad(20.0), timetorotang(time_s, 1.))
    end
    partial_pix_idx = Healpix.ang2pixRing.(Ref(Healpix.Resolution(NSIDE)), dirs[:, 1], dirs[:, 2])
    global pix_idx = append!(pix_idx, partial_pix_idx)
end

noise_tod = generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, total_time, fsamp_hz, slope, sigma_k, fknee_hz, rank, comm, 1234)
(destr_map, a) = destripe(pix_idx, noise_tod, num_of_pixels, data_properties, rank, comm; threshold = 1e-9, max_iter = 1000)

# 1 polarimeter
num_of_polarimeters = 1
baselines_per_process = Int64(total_time / baseline_length_s) * num_of_polarimeters
sigma_k =  [0.001]
fknee_hz = [0.01]
slope = [-1]

chunks = Any[Any[datachunk(1, 1, 6, 6)]]
data_properties = [
    TodNoiseProperties(pol = 1, rms = 0.001, baselines = [100, 100, 100, 100, 100, 100]),
]

pix_idx = Int64[]

for i in 1:1  #loop on detectors
    times = 0.0:tau_s:total_time - tau_s

    (dirs, ψ) = genpointings([0.,0.,1.], times; latitude_deg = 28.29) do time_s
        return (0.0, deg2rad(20.0), timetorotang(time_s, 1.))
    end

    partial_pix_idx = Healpix.ang2pixRing.(Ref(Healpix.Resolution(NSIDE)), dirs[:, 1], dirs[:, 2])
    global pix_idx = append!(pix_idx, partial_pix_idx)

end

noise_tod = generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, total_time, fsamp_hz, slope, sigma_k, fknee_hz, rank, comm, 1234)

(destr_map_sigma1, a_1) = destripe(pix_idx, noise_tod, num_of_pixels, data_properties, rank, comm; threshold = 1e-9, max_iter = 1000)

diff = destr_map[isfinite.(destr_map)] - destr_map_sigma1[isfinite.(destr_map_sigma1)]
@test sum(diff .- mean(diff)) ≈ 0. atol = 1e-10