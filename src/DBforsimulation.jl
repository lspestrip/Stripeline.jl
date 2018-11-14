export get_info_from_DB
import Random


"""
    function get_info_from_DB(db, horns)
    
        This function can be used to take information about specific horns from the STRIP instrument database.

        Given:
        -The database (db = Stripeline.InstrumentDB())
        -A vector containg the names of the horns of interest  (e.g. ["Y0", "I1", "R6"])

        It returns 5 arrays, containing for each horn:
        - the 3D vector containing the orientation of the horn in the sky
        - the unique ID of the polarimeter associated to the horn
        - the estimated bandwidth, in Hz
        - the estimated noise temperature, in K
        - the estimated knee frequency, in Hz
        
        For those polarimeters without an available measurement of the knee frequency (13 out of 49)
        the function randomly picks a value from a vector containing all the 36 fknee measurements available in the database.

"""
function get_info_from_DB(db, horns)
    orientations=[]
    polarimeters = Array{Int64}(undef, length(horns))
    β_hz = Array{Float64}(undef,length(horns))
    tnoise_k = Array{Float64}(undef,length(horns))
    fknee_hz = Array{Float64}(undef,length(horns))

    measured_fknees = [0.13, 0.04, 1.3, 0.06, 0.23, 0.025, 0.016,
    0.9, 0.12, 0.3, 0.2, 0.09, 0.31, 0.05, 
    0.1, 0.02, 0.027, 0.07, 0.02, 0.014, 0.007,
    0.41, 0.4, 0.09, 0.138, 0.05, 0.17, 0.6, 
    0.028, 0.0048, 0.01, 0.058, 0.01, 0.08, 
    0.016, 0.09]

    rng = Random.MersenneTwister(1234)

    for i in 1:length(horns)
        append!(orientations, [db.focalplane[horns[i]].orientation])
        polarimeters[i] = db.focalplane[horns[i]].polid
        β_hz[i] = db.detectors[polarimeters[i]].bandshape.bandwidth_hz
        tnoise_k[i] = db.detectors[polarimeters[i]].tnoise.tnoise_k
          
        if(db.detectors[polarimeters[i]].spectrum.fknee_q_hz ==0)
            fknee_hz[i] = rand(rng, measured_fknees)
        else
            fknee_hz[i] = db.detectors[polarimeters[i]].spectrum.fknee_q_hz
        end
    end

    return orientations, polarimeters, β_hz, tnoise_k, fknee_hz

end


