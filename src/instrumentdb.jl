export Horn, Detector, InstrumentDB, BandshapeInfo, SpectrumInfo, NoiseTemperatureInfo

import YAML

doc"""
Information about a STRIP horn

This structure holds a number of parameters relative to each feed horn in the
STRIP focal plane.

You should initialize `Horn` objects via the `InstrumentDB` constructor, which
loads their definition from a STRIP instrument database in YAML format.

Name                 | Type           | Meaning
-------------------- | -------------- | -------------------------------------------------------
`name`               | String         | Name of the horn, e.g., `I0`
`id`                 | Int            | Unique number of the horn, starting from 1
`polid`              | Int            | Unique ID of the polarimeter associated with the horn
`moduleid`           | Int            | Number of the horn within the module, from 0 to 6
`color`              | String         | Name of the color associated with the module
`orientation`        | Array{Float64} | 3D vector containing the orientation of the horn in the sky
`fwhm_x_deg`         | Float64        | FWHM of the beam along the X axis, in degrees
`fwhm_y_deg`         | Float64        | FWHM of the beam along the Y axis, in degrees
`main_spillover`     | Float64        | Main reflector spillover
`sub_spillover`      | Float64        | Sub-reflector spillover
`xpd_db`             | Float64        | Cross-polarization, in dB
`directivity_dbi`    | Float64        | Directivity, in dBi
`ellipticity`        | Float64        | Ellipticity
"""
struct Horn
    name::String
    id::Int
    polid::Int
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

doc"""
    BandshapeInfo

Information about the spectral band response of a polarimeter.
"""
struct BandshapeInfo
    center_frequency_hz::Float64
    center_frequency_err_hz::Float64
    bandwidth_hz::Float64
    bandwidth_err_hz::Float64
    lowest_frequency_hz::Float64
    highest_frequency_hz::Float64
    num_of_frequencies::Int
    response::Array{Float64,1}
    test_id::Int
    analysis_id::Int
end

doc"""
    BandshapeInfo()

Initialize a BandshapeInfo object with all values set to zero.
"""
BandshapeInfo() = BandshapeInfo(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, Float64[], 0, 0)

doc"""
    SpectrumInfo

Information about the noise spectrum of the output of a polarimeter.
"""
struct SpectrumInfo
    slope_i::Float64
    slope_q::Float64
    slope_u::Float64
    slope_err_i::Float64
    slope_err_q::Float64
    slope_err_u::Float64
    fknee_i_hz::Float64
    fknee_q_hz::Float64
    fknee_u_hz::Float64
    fknee_i_err_hz::Float64
    fknee_q_err_hz::Float64
    fknee_u_err_hz::Float64
    wn_i_k2_hz::Float64
    wn_q_k2_hz::Float64
    wn_u_k2_hz::Float64
    wn_i_err_k2_hz::Float64
    wn_q_err_k2_hz::Float64
    wn_u_err_k2_hz::Float64
    test_id::Int
    analysis_id::Int
end

doc"""
    SpectrumInfo()

Initialize a SpectrumInfo object with all values set to zero.
"""
SpectrumInfo() = SpectrumInfo(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)

doc"""
    NoiseTemperatureInfo

Information about the noise temperature of a polarimeter. This structure is used
for the field `tnoise` of the `Detector` struct. """
struct NoiseTemperatureInfo
    tnoise_k::Float64
    tnoise_err_k::Float64
    tnoise_test_ids::Array{Int,1}
    tnoise_analysis_ids::Array{Int,1}
    tnoise_values_k::Array{Float64,1}
end

doc"""
    NoiseTemperatureInfo()

Initialize a NoiseTemperatureInfo object with all values set to zero.
"""
NoiseTemperatureInfo() = NoiseTemperatureInfo(0.0, 0.0, Int[], Int[], Float64[])

doc"""
Information about a STRIP detector

This structure holds information about a STRIP polarimeter.

You should initialize `Detector` objects via the `InstrumentDB` constructor,
which loads their definition from a local STRIP instrument database.
"""
struct Detector
    id::Int
    name::String
    band::String
    
    bandshape::BandshapeInfo
    spectrum::SpectrumInfo
    tnoise::NoiseTemperatureInfo
end

doc"""
STRIP instrument database

The "database" contains information about feed horns (through the structure
`Horn`) and polarimeters (through the structure `Detector`). You should usually
create an object of this kind using the default constructor, which parses a
set of YAML files containing the real parameters of the instrument.
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

function parsebandshape(banddict::Dict{Any,Any})
    BandshapeInfo(get(banddict, "center_frequency_hz", 0.0),
        get(banddict, "center_frequency_err_hz", 0.0),
        get(banddict, "bandwidth_hz", 0.0),
        get(banddict, "bandwidth_err_hz", 0.0),
        get(banddict, "lowest_frequency_hz", 0.0),
        get(banddict, "highest_frequency_hz", 0.0),
        get(banddict, "num_of_frequencies", 0),
        get(banddict, "response", Float64[]),
        get(banddict, "test_id", 0),
        get(banddict, "analysis_id", 0))
end

function parsespectrum(specdict::Dict{Any,Any})
    SpectrumInfo(get(specdict, "I_slope", 0.0),
        get(specdict, "Q_slope", 0.0),
        get(specdict, "U_slope", 0.0),
        get(specdict, "I_slope_err", 0.0),
        get(specdict, "Q_slope_err", 0.0),
        get(specdict, "U_slope_err", 0.0),
        get(specdict, "I_fknee_hz", 0.0),
        get(specdict, "Q_fknee_hz", 0.0),
        get(specdict, "U_fknee_hz", 0.0),
        get(specdict, "I_fknee_err_hz", 0.0),
        get(specdict, "Q_fknee_err_hz", 0.0),
        get(specdict, "U_fknee_err_hz", 0.0),
        get(specdict, "I_wn_level_k2_hz", 0.0),
        get(specdict, "Q_wn_level_k2_hz", 0.0),
        get(specdict, "U_wn_level_k2_hz", 0.0),
        get(specdict, "I_wn_level_err_k2_hz", 0.0),
        get(specdict, "Q_wn_level_err_k2_hz", 0.0),
        get(specdict, "U_wn_level_err_k2_hz", 0.0),
        get(specdict, "test_id", 0),
        get(specdict, "analysis_id", 0))
end

function parsetnoise(tnoisedict::Dict{Any,Any})
    NoiseTemperatureInfo(get(tnoisedict, "tnoise_k", 0.0),
        get(tnoisedict, "tnoise_err_k", 0.0),
        get(tnoisedict, "tnoise_test_ids", Int[]),
        get(tnoisedict, "analysis_ids", Int[]),
        get(tnoisedict, "values_k", Float64[]))
end

function parsedetdict(detdict)
    detectors = Dict{Int,Detector}()
    for curdet in detdict
        id = curdet["id"]

        bandpass = haskey(curdet, "bandpass") ? parsebandshape(curdet["bandpass"]) : BandshapeInfo()
        spectrum = haskey(curdet, "spectrum") ? parsespectrum(curdet["spectrum"]) : SpectrumInfo()
        tnoise = haskey(curdet, "tnoise") ? parsetnoise(curdet["tnoise"]) : NoiseTemperatureInfo()

        detectors[id] = Detector(id,
                                 get(curdet, "name", @sprintf("STRIP%02d", id)),
                                 get(curdet, "band", id ≤ 70 ? "Q" : "W"),
                                 bandpass,
                                 spectrum,
                                 tnoise)
    end
    
    detectors
end

doc"""
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

doc"""
    InstrumentDB()

Load the default STRIP instrument database. Return an instance of
a InstrumentDB object."""
InstrumentDB() = InstrumentDB(joinpath(Pkg.dir("Stripeline"), "instrumentdb"))
