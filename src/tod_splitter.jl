export datachunk
export split_into_n, split_tod_mpi, get_chunk_properties, generate_noise_mpi

import Healpix
import CorrNoise
using Random
using FITSIO


try
    import MPI
catch
end 


"""
This structure holds a number of parameters relative to a certain chunk of data.

Field              | Type           | Meaning
:----------------- |:-------------- |:------------------------------------------------------------
`pol_number`       | Int            | ID number of the polarimeter
`first_idx`        | Int            | Index of the first element of the chunk
`last_idx`         | Int            | Index of the last element of the chunk
`num_of_elements`  | Int            | Total number of elements of the chunk

# Example
given the following array of data: [1, 10, 20, 3, 4, 5, 7]
measured by the polarimeter number 2

the chunk: [20, 3, 4]
will correspond to the following structure:
datachunk(2, 3, 5, 3)
"""
struct datachunk
    pol_number ::Int
    first_idx::Int
    last_idx::Int
    num_of_elements::Int
end 


"""
    function split_into_n(length, num_of_segments)
        Given the `length` of the array, it convenientely 
        splits it into `num_of_segments` sections of as similar length as possible.
        It returns an array containing the number of elements of each section.

        # Example
        julia> split_into_n(20, 3)
        3-element Array{Int64,1}:
        6
        7
        7
"""
function split_into_n(length, num_of_segments)
    @assert num_of_segments >0
    @assert length >= num_of_segments
    start_pos = zeros(Int, num_of_segments+1)
    
    for i in 1:num_of_segments+1
        start_pos[i] =  floor(((i-1)*length/num_of_segments))
    end
   
    return start_pos[2:end]-start_pos[1:end-1]    
end


"""
    function split_tod_mpi(total_time, baseline_length_s, baselines_per_process, num_of_MPI_proc)
           
        This function can be used to split the TOD production of many polarimeters among MPI processes.

        It requires in input:

        -the total time (in seconds) of the simulated observation 
        -the length (in seconds) of each 1/f noise baseline
        -the array containing the number of 1/f baselines to simulate for each process.
            It can be obtained by using the function `plit_into_n` in the following way:

            split_into_n(num_of_polarimeters*baselines_per_pol, num_of_MPI_proc)
                
            where baselines_per_pol = Int64(total_time/baseline_length_s) is the number of baselines of each polarimeter

        -the number of MPI processes used

        It returns an array of arrays of `datachunk` instances, of length == num_of_MPI_proc
        where each element tells the chunk of data that each process should simulate (see Example) 

        # Example
        ```julia-repl
        julia> num_of_polarimeters = 4
        julia> num_of_MPI_proc = 3
        julia> total_time = 50
        julia> baseline_length_s = 10
        julia> baselines_per_pol = Int64(total_time/baseline_length_s)
        5

        julia> baselines_per_process = split_into_n(20, 3)
        3-element Array{Int64,1}:
        6
        7
        7

        julia> chunks = split_tod_mpi(total_time, baseline_length_s, baselines_per_process, num_of_MPI_proc)
        3-element Array{Any,1}:
        Any[datachunk(1, 1, 5, 5), datachunk(2, 1, 1, 1)]
        Any[datachunk(2, 2, 5, 4), datachunk(3, 1, 3, 3)]
        Any[datachunk(3, 4, 5, 2), datachunk(4, 1, 5, 5)]
        ```

        which means:

        - process number 0 should simulate: 
        polarimeter number 1 from baseline 1 to baseline 5,  total number of baselines = 5
        polarimeter number 2 from baseline 1 to baseline 1 , total number of baselines = 1

        - process number 1 should simulate: 
        polarimeter number 2 from baseline 2 to baseline 5,  total number of baselines = 4
        polarimeter number 3 from baseline 1 to baseline 3 , total number of baselines = 3

        - process number 2 should simulate: 
        polarimeter number 3 from baseline 4 to baseline 5,  total number of baselines = 2
        polarimeter number 4 from baseline 1 to baseline 5,  total number of baselines = 5

"""
function split_tod_mpi(total_time, baseline_length_s, baselines_per_process, num_of_MPI_proc)

    duration = Int64(total_time/baseline_length_s)

    #initialization
    detector_num = 1
    sample_idx = 0
    samples_in_det = duration
    result = []
    
    for rank_num in 0:(num_of_MPI_proc-1)  #loop on MPI processes
        
        samples_for_this_process = baselines_per_process[rank_num+1]
        samples_left = samples_for_this_process
        data_this_rank = []
        
        while samples_left > 0  #loop on detectors
        
            #if the current detector has more samples than needed to fill the current MPI process
            if samples_in_det > samples_left
                
                first_idx = sample_idx+1
                last_idx = sample_idx+samples_left
                data = datachunk(detector_num, first_idx, last_idx, samples_left)
                data_this_rank = append!(data_this_rank, [data])
                
                sample_idx = sample_idx + samples_left
                samples_in_det = samples_in_det - samples_left
                samples_left = 0 
            
            #if the current detector has not enough samples to provide the current MPI process 
            #with the required number of samples. In this case we need to increase "detector_num" before the next iteration
            else 
                
                first_idx = sample_idx+1
                last_idx = sample_idx+samples_in_det
                data = datachunk(detector_num, first_idx, last_idx, samples_in_det)
                data_this_rank = append!(data_this_rank, [data])
               
                samples_left = samples_left - samples_in_det
                samples_in_det = 0
            end
            
            if samples_in_det == 0
                detector_num +=1
                sample_idx = 0
                samples_in_det = duration 
            end
        
        end   
        
        result = append!(result, [data_this_rank])
        
    end 
        
    return result
    
end

"""
    function get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)
    
        Given:
        - the data chunks (which can be obtained by using the function `split_tod_mpi`)
        - the length (in seconds) of each 1/f noise baseline
        - the sampling frequency (in Hz)
        - the number of current MPI rank

        this function extracts useful information to perform the TOD simulation in the current rank.


        It returns a tuple containing 5 arrays:
        - the ID number of the polarimeters that the current rank will simulate
        - the start time of the acquisition portion for each polarimeter
        - the stop time of the acquisition portion for each polarimeter
        - the number of 1/f baselines for each polarimeter
        - the total number of samples for each polarimeter

        # Example
        ```julia-repl
        julia> baseline_length_s = 10
        julia> fsamp_hz = 10
        julia> rank = 1

        julia> chunks = [[datachunk(1, 1, 5, 5), datachunk(2, 1, 1, 1)], [datachunk(2, 2, 5, 4), datachunk(3, 1, 3, 3)], [datachunk(3, 4, 5, 2), datachunk(4, 1, 5, 5)]]

        get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)
        ([2, 3], [10.0, 0.0], [50.0, 30.0], [4, 3], [400, 300])

        which means that rank 1 will simulate:
        - polarimeter number 2 from 10 s (from the starting of the acquisition) to 50 s, 
          with a total of 4 1/f baselines and 400 samples.
        - polarimeter number 3 from 0 s (from the starting of the acquisition) to 30 s, 
          with a total of 3 1/f baselines and 300 samples.

    """
    function get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)

        this_rank_chunk = chunks[rank+1]
        first_time, last_time = [Array{Float64}(undef, length(this_rank_chunk)) for i in (1:2)]
        detector_number, num_of_baselines, baseline_len, num_of_samples = [Array{Int64}(undef, length(this_rank_chunk)) for i in (1:4)]
    
        for i in 1:length(this_rank_chunk)
            detector_number[i] = this_rank_chunk[i].pol_number
            first_time[i] = (this_rank_chunk[i].first_idx-1)*baseline_length_s
            last_time[i] = this_rank_chunk[i].last_idx*baseline_length_s -0.99*(1/fsamp_hz)
            num_of_baselines[i] =  this_rank_chunk[i].num_of_elements
            num_of_samples[i] = num_of_baselines[i]*baseline_length_s*fsamp_hz  
        end
        return (detector_number, first_time, last_time, num_of_baselines, num_of_samples)
    end
    

"""
    function generate_noise_mpi(chunks, baseline_length_s, baselines_per_process, fsamp_hz, σ_k, fknee_hz, rank, comm)

    This MPI based function is useful to generate white noise and 1/f noise when the simulation is splitted in various MPI processes.

    It requires in input:
    - the data chunks (which can be obtained by using the function `split_tod_mpi`)
    - the length (in seconds) of each 1/f noise baseline
    - the array containing the number of 1/f baselines to simulate for each process.
            It can be obtained by using the function `plit_into_n` in the following way:

            split_into_n(num_of_polarimeters*baselines_per_pol, num_of_MPI_proc)
                
            where baselines_per_pol = Int64(total_time/baseline_length_s) is the number of baselines of each polarimeter

    - the sampling frequency (in Hz)
    - The RMS of the temperature (in K), σ_k = tsys_k / sqrt(β_hz * τ_s)
    - The knee frequency (in Hz)
    - the MPI communicator

    It returns the noise TOD for the current rank. 

    Each partial noise TOD is generated by rank 0 and then sent to the correspondent rank. 

    N.B. if you want to use this function without MPI, remember to put rank = 0 and comm = Nullable{}()


"""
function generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, fsamp_hz, σ_k, fknee_hz, rank, comm)

    if rank == 0
        previous_detector = chunks[1][1].pol_number
        noise = Float64[]  #noise rank 0
        seed = rand(1:1000)
        rng = CorrNoise.OofRNG(CorrNoise.GaussRNG(MersenneTwister(seed)), -1, 1.15e-5, fknee_hz, fsamp_hz)
    
        for i in 1:length(chunks)   #loop on ranks
            num_noise_samples = baselines_per_process[i]*baseline_length_s*fsamp_hz

            for j in 1:length(chunks[i]) #loop on detectors per rank
                cur_detector = chunks[i][j].pol_number

                #if cur_detector != previous_detector  #if new detector generate new noise
                    #rng = CorrNoise.OofRNG(CorrNoise.GaussRNG(MersenneTwister(1234)), -1, 1.15e-5, fknee_hz, fsamp_hz)
                #end

                cur_num_noise_samples = chunks[i][j].num_of_elements*baseline_length_s*fsamp_hz
                cur_noise = Float64[CorrNoise.randoof(rng) * σ_k for i in 1:(cur_num_noise_samples)]
                previous_detector =  cur_detector

                if (i > 1)  #send to other ranks, not rank 0
                        MPI.Send(cur_noise, i-1, 0, comm)
                else
                    noise = append!(noise, cur_noise)
                end
            end   
        end
    else
        noise = Float64[]
        for j in 1:length(chunks[rank+1])
            cur_noise = Array{Float64}(undef, chunks[rank+1][j].num_of_elements*baseline_length_s*fsamp_hz)
            MPI.Recv!(cur_noise, 0, 0, comm)
            noise = append!(noise, cur_noise)
        end
    end

    return noise
end

