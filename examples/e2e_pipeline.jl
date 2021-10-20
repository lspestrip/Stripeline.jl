import MPI
import Healpix
import Random
import CorrNoise
import Stripeline
const Sl = Stripeline

using Printf

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
commsize = MPI.Comm_size(comm)

printmsg(msg) = (rank == 0) && print(msg)

printmsg(@sprintf("""MPI parameters:

- Rank: %d
- Number of MPI processes: %d

""", rank, commsize))

days_of_observation = 1
num_of_polarimeters = 1
baseline_length_s = 10
fsamp_hz = 100
NSIDE = 256

printmsg(@sprintf(
    """Parameters of the simulation:

- Days of observation: %d
- Number of polarimeters: %d
- Length of each baseline: %.1f s
- Sampling frequency: %.1f Hz

""",
    days_of_observation,
    num_of_polarimeters,
    baseline_length_s,
    fsamp_hz,
))

tcmb_k = 3
tatm_k = 15
ttel_k = 3
tnoise_k = 35
fknee_hz = 0.01
β_hz = 7e9

total_time_s = days_of_observation * 24 * 3600.0
tsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k
τ_s = 1 / fsamp_hz
σ_k = (tsys_k / sqrt(β_hz * τ_s))

samples_per_pol = Int64(total_time_s * fsamp_hz)
baselines_per_pol = Int(total_time_s / baseline_length_s)

sky_map = joinpath("test", "testfiles", "PySM_inputmap_nside256.fits")

printmsg("Reading map \"$(sky_map)\"\n")
inputmap = Healpix.readMapFromFITS(
    sky_map,
    2,
    Float32,
)
inputmap_resol = inputmap.resolution

resol = Healpix.Resolution(NSIDE)
num_of_pixels = resol.numOfPixels

baselines_per_process = Sl.split_into_n(num_of_polarimeters * baselines_per_pol, commsize)
chunks = Sl.split_tod_mpi(total_time_s, baseline_length_s, baselines_per_process, commsize)
this_rank_chunk = chunks[rank + 1]

(detector_number, first_time, last_time, num_of_baselines, num_of_samples) = Sl.get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)

printmsg("Projecting the map onto the TOD\n")

pix_idx = Int64[]
pix_idx_inputmap = Int64[]
sky_tod = Float32[]

for i in 1:length(this_rank_chunk)
    local times = first_time[i]:τ_s:last_time[i]

    local (dirs, ψ) = Sl.genpointings(
        Float32[0.0, 0.0, 1.0],
        times;
        latitude_deg = 28.29,
    ) do time_s
        (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1.))
    end

    partial_pix_idx_inputmap = Healpix.ang2pixRing.(
        Ref(inputmap_resol),
        dirs[:, 1],
        dirs[:, 2],
    )
    partial_pix_idx = Healpix.ang2pixRing.(
        Ref(resol),
        dirs[:, 1],
        dirs[:, 2],
    )

    partial_sky_tod = inputmap.pixels[partial_pix_idx_inputmap]
    global sky_tod = append!(sky_tod, partial_sky_tod)
    global pix_idx = append!(pix_idx, partial_pix_idx)
end

printmsg("Generating the noise\n")

noise_tod = Sl.generate_noise_mpi(
    chunks,
    baselines_per_process,
    baseline_length_s,
    total_time_s,
    fsamp_hz,
    σ_k,
    fknee_hz,
    1.0,
    rank = rank,
    comm = comm,
)

tod = sky_tod + noise_tod

printmsg("Running the destriper\n")

function callback(
    ;
    iter_idx,
    max_iter,
    convergence_parameter,
    convergence_threshold,
)
    @printf("%d/%d, %e > %e\r", iter_idx, max_iter, convergence_parameter, convergence_threshold)
end
println()

baseline_len = repeat([baseline_length_s * fsamp_hz], baselines_per_process[rank + 1])
data_properties = Sl.build_noise_properties([1], [σ_k], num_of_baselines, num_of_samples)
results = Sl.destripe(pix_idx, tod, num_of_pixels, data_properties, rank, comm = comm, callback = callback)

out_map_name = "destriped_map.fits"
printmsg("Saving the map in \"$(out_map_name)\"\n")

if rank == 0
    mapfile = Healpix.HealpixMap{Float64, Healpix.RingOrder}(NSIDE)
    mapfile.pixels = results.best_sky_map

    Healpix.saveToFITS(mapfile, out_map_name, typechar = "D")
end

MPI.Finalize()
