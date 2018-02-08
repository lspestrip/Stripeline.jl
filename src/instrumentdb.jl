export Horn, Detector, InstrumentDB

import YAML

"""
    Information about a STRIP horn

You should initialize Horn objects via the InstrumentDB constructor,
which loads their definition from a local STRIP instrument database.
"""
struct Horn
    name::String
    id::Int
    polid::String
    moduleid::Int
    color::String
    orientation::Array{Float64,1}
    fwhm_x_deg::Float64
    fwhm_y_deg::Float64
    main_spillover::Float64
    sub_spillover::Float64
    xpd_db::Float64
    directivity_dbi::Float64
    ellipticity::Float64
end

"""
    Information about a STRIP detector

You should initialize Detector objects via the InstrumentDB constructor,
which loads their definition from a local STRIP instrument database.
"""
struct Detector
    id::Int
    name::String
    center_frequency_hz::Float64
    bandwidth_hz::Float64
    alpha_i::Float64
    alpha_q::Float64
    alpha_u::Float64
    fknee_q_hz::Float64
    fknee_u_hz::Float64
    wn_q_k2_hz::Float64
    wn_u_k2_hz::Float64
end

"""
    STRIP instrument database
"""
struct InstrumentDB
    focalplane::Dict{String,Horn}
    detectors::Dict{Int,Detector}
end

function parsefpdict(fpdict::Dict{Any,Any})
    focalplane = Dict{String,Horn}()
    for (key, value) in fpdict["horns"]
        refhorn = if haskey(fpdict["pairs"], key)
            refhornname = fpdict["pairs"][key]
            fpdict["horns"][refhornname]
        else
            value
        end
        
        focalplane[key] = Horn(key,
                               value["id"],
                               value["polarimeter_id"],
                               value["module_id"],
                               value["color"],
                               value["orientation"],
                               refhorn["fwhm_x_deg"],
                               refhorn["fwhm_y_deg"],
                               refhorn["main_spillover"],
                               refhorn["sub_spillover"],
                               refhorn["xpd_db"],
                               refhorn["directivity_dbi"],
                               refhorn["ellipticity"])
    end
    
    focalplane
end

function parsedetdict(detdict::Dict{Any,Any})
    detectors = Dict{Int,Detector}()
    for curdet in detdict["polarimeters"]
        id = curdet["id"]
        detectors[id] = Detector(id,
                                 curdet["name"],
                                 get(curdet, "center_frequency_hz", 0.0),
                                 get(curdet, "bandwidth_hz", 0.0),
                                 get(curdet, "alpha_i", 0.0),
                                 get(curdet, "alpha_q", 0.0),
                                 get(curdet, "alpha_u", 0.0),
                                 get(curdet, "fknee_q_hz", 0.0),
                                 get(curdet, "fknee_u_hz", 0.0),
                                 get(curdet, "wn_q_k2_hz", 0.0),
                                 get(curdet, "wn_q_k2_hz", 0.0))
    end
    
    detectors
end

"""
    InstrumentDB(dbpath)

Load the STRIP instrument database from the specified path.
Return an instance of a InstrumentDB object.
"""
function InstrumentDB(dbpath)
    focalplanedict = open(joinpath(dbpath, "strip_focal_plane.yaml")) do f
        YAML.load(f)
    end
    
    detectordict = open(joinpath(dbpath, "strip_detectors.yaml")) do f
        YAML.load(f)
    end
    
    InstrumentDB(parsefpdict(focalplanedict), parsedetdict(detectordict))
end
