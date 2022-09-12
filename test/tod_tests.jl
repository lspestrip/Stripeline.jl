using AstroTime
import Stripeline
const Sl = Stripeline
import Statistics: cov

function test_split_into_n()
    @test split_into_n(20, 3)  == [6, 7, 7]
end

function test_allocate_tod(time_range)
    t0 = minimum(time_range)
    t1 = maximum(time_range)

    # Number of samples acquired during the observation (this is usually huge!)
    nsamples = length(time_range)

    # Three fake polarimeters
    polarimeters = 1:3

    mpi_size = 10  # Number of fake MPI processes

    # Allocate all the TODs in memory, assuming we are running `mpi_size` processes
    tods = [allocate_tod(
        StripTod,
        Float32,
        time_range,
        polarimeters,
        mpi_rank = rank,
        mpi_size = mpi_size,
        rng_seed = 12345,
    ) for rank in 0:(mpi_size - 1)]

    # Have the N TODs been really allocated?
    @test length(tods) == 10

    # The first TOD must start when the observation begins…
    @test minimum(tods[1].time_range) == t0
    # …and the last TOD must end when the observation ends
    @test maximum(tods[end].time_range) == t1

    # Check that the total number of samples in each TOD matches what we have asked for
    @test sum([length(cur_tod.time_range) for cur_tod in tods]) == length(time_range)

    # Check that there are no overlaps between consecutive TODs
    for i in 2:length(tods)
        @test maximum(tods[i - 1].time_range) < minimum(tods[i].time_range)
    end

    for i in eachindex(tods)
        cur_shape = size(tods[i].samples)
        @test cur_shape[1] == length(tods[i].time_range)
        @test cur_shape[2] == 8
        @test cur_shape[3] == length(polarimeters)
    end

    # Check that the random number generator generates a different
    # sequence across the MPI processes
    @test rand(tods[1].rng, Int8) == 67
    @test rand(tods[2].rng, Int8) == -109
    @test rand(tods[3].rng, Int8) == 47
end


function test_stokes()
    # Three fake polarimeters
    polarimeters = 1:3

    mpi_size = 10  # Number of fake MPI processes

    # Allocate all the TODs in memory, assuming we are running `mpi_size` processes
    tods = [allocate_tod(
        StripTod,
        Float32,
        0.0:0.1:10.0,
        polarimeters,
        mpi_rank = rank,
        mpi_size = mpi_size,
        rng_seed = 12345,
    ) for rank in 0:(mpi_size - 1)]

    stokestods = [stokes(tod) for tod in tods]

    for (tod, stokestod) in zip(tods, stokestods)
        @assert size(stokestod.samples)[1] == size(tod.samples)[1]
        @assert size(stokestod.samples)[2] == 3
        @assert size(stokestod.samples)[3] == size(tod.samples)[3]
    end
end


function test_fillnoise()

    db = Sl.InstrumentDB()
    tod = Sl.allocate_tod(StripTod, Float32, 0.0:0.1:10000.0, ["I3"])
    Sl.fillnoise!(tod) do polname
        s = Sl.spectrum(db, polname)
        (s.pwr_cov_matrix_k2, s.dem_cov_matrix_k2)
    end

    expected = Sl.spectrum(db, tod.polarimeters[1]).pwr_cov_matrix_k2
    computed = cov(tod.samples[:, Sl.PWR_Q1_RANK:Sl.PWR_U2_RANK, 1])

    # PWR series have a higher covariance
    @test all(@. expected - computed < 1e-1)

    expected = Sl.spectrum(db, tod.polarimeters[1]).dem_cov_matrix_k2
    computed = cov(tod.samples[:, Sl.DEM_Q1_RANK:Sl.DEM_U2_RANK, 1])

    # DEM series have a lower covariance
    @test all(@. expected - computed < 1e-5)

end

test_split_into_n()

# Run tests using a simple floating-point range
test_allocate_tod(0.0:0.01:1000.0)

# Re-run `test_allocate_tod` using an AstroTime range
let fsamp_hz = 100.0,
    τ_s = 1 / fsamp_hz,                       # Integration time for a sample
    t0 = TAIEpoch(2022, 1, 1, 12, 0, 0.0),    # Fake date when the observation starts
    t1 = t0 + 1hours,                         # End of the observation
    time_range = t0 : (τ_s * seconds) : t1    # Range covering all the samples

    test_allocate_tod(time_range)
end

test_stokes()

test_fillnoise()
