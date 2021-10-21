using AstroTime
import Stripeline
const Sl = Stripeline

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
        Float32,
        time_range,
        polarimeters,
        mpi_rank = rank,
        mpi_size = mpi_size,
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
end

test_split_into_n()

# Run tests using a simple floating-point range
test_allocate_tod(0.0:0.01:1000.0)

# Re-run `test_allocate_tod` using an AstroTime range
begin
    fsamp_hz = 100.0
    τ_s = 1 / fsamp_hz  # Integration time for a sample

    # Fake date when the observation starts
    t0 = TAIEpoch(2022, 1, 1, 12, 0, 0.0)
    # End of the observation
    t1 = t0 + 1hours

    # This is a range that covers *all* the samples during the observation
    time_range = t0 : (τ_s * seconds) : t1
    
    test_allocate_tod(time_range)
end

