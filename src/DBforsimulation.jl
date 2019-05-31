export get_info_from_DB
import Random


"""
    function get_info_from_DB(db, horns)
    
Given an instance of a [`InstrumentDB`](@ref) object and a vector containing
horn names (e.g., `["Y0", "I1", "R6"]`), this function returns a tuple
containing the following fields:

1. An array of 3D vectors containing the orientation of the horn in the sky
2. An array of integers containing the unique ID of each polarimeter associated
   to the horn
3. An array of floats containing the estimated bandwidth, in Hz
4. An array of floats containing the estimated noise temperature, in K
5. An array of floats containing the estimated knee frequency, in Hz
        
For those polarimeters without an available measurement of the knee frequency
(13 out of 49) the function randomly picks a value from a vector containing all
the 36 fknee measurements available in the database.
"""
function get_info_from_DB(db, horns)
    orientations = []
    polarimeters = Array{Int64}(undef, length(horns))
    β_hz = Array{Float64}(undef, length(horns))
    tnoise_k = Array{Float64}(undef, length(horns))
    fknee_hz = Array{Float64}(undef, length(horns))

    measured_fknees = Float64[]
    for horn in keys(db.focalplane)
        det = db.detectors[db.focalplane[horn].polid]
        det.spectrum.fknee_q_hz > 0 && push!(measured_fknees, det.spectrum.fknee_q_hz)
        det.spectrum.fknee_u_hz > 0 && push!(measured_fknees, det.spectrum.fknee_u_hz)
    end

    rng = Random.MersenneTwister(1234)

    for i in 1:length(horns)
        append!(orientations, [db.focalplane[horns[i]].orientation])
        polarimeters[i] = db.focalplane[horns[i]].polid
        β_hz[i] = db.detectors[polarimeters[i]].bandshape.bandwidth_hz
        tnoise_k[i] = db.detectors[polarimeters[i]].tnoise.tnoise_k
          
        if(db.detectors[polarimeters[i]].spectrum.fknee_q_hz > 0)
            fknee_hz[i] = db.detectors[polarimeters[i]].spectrum.fknee_q_hz
        else
            fknee_hz[i] = rand(rng, measured_fknees)
        end
    end

    (orientations, polarimeters, β_hz, tnoise_k, fknee_hz)
end
