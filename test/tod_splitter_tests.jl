using Test
using Statistics

num_of_polarimeters = 4
num_of_MPI_proc = 3
total_time = 50
baseline_length_s = 10
fsamp_hz = 10
fknee_hz = [0.01, 0.01]
σ_k = [0.1, 100]
rank = 1

baselines_per_process = [6, 7, 7]

@test split_tod_mpi(
    total_time,
    baseline_length_s,
    baselines_per_process,
    num_of_MPI_proc,
) == [
    [datachunk(1, 1, 5, 5), datachunk(2, 1, 1, 1)],
    [datachunk(2, 2, 5, 4), datachunk(3, 1, 3, 3)],
    [datachunk(3, 4, 5, 2), datachunk(4, 1, 5, 5)],
]

chunks = [
    [datachunk(1, 1, 5, 5), datachunk(2, 1, 1, 1)],
    [datachunk(2, 2, 5, 4), datachunk(3, 1, 3, 3)],
    [datachunk(3, 4, 5, 2), datachunk(4, 1, 5, 5)],
]

@test get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank) ==
      ([2, 3], [10.0, 0.0], [49.901, 29.901], [4, 3], [400, 300])

