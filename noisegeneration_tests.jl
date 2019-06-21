using Test
using Statistics


baseline_length_s = 10
fsamp_hz = 10
fknee_hz = [0.01, 0.01]
σ_k = [0.1, 100]
slope = [-1,-1]


chunks = [[datachunk(1, 1, 1000, 1000), datachunk(2, 1, 1000, 1000)]]
baselines_per_process = 2000
rank=0
comm = missing

noise =  generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, total_time, fsamp_hz, σ_k, fknee_hz, slope,rank, comm, 1234)
@test length(noise) == fsamp_hz*baseline_length_s*baselines_per_process
@test std(noise[1:100000]) ≈ σ_k[1]  rtol=0.01
@test std(noise[end-100000: end]) ≈ σ_k[2]  rtol=0.01