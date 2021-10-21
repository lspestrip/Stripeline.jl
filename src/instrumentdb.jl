export Horn, Detector, InstrumentDB, BandshapeInfo, SpectrumInfo, NoiseTemperatureInfo
export InstrumentDB, defaultdbfolder, parsefpdict, parsedetdict
export sensitivity_tant, t_to_trj, trj_to_t, deltat_to_deltatrj, deltatrj_to_deltat
export detector, bandpass, bandshape, spectrum, fknee_hz, tnoise

import YAML
import Stripeline
import Base: show
import RecipesBase

using Printf

@doc raw"""
Information about a STRIP horn

This structure holds a number of parameters relative to each feed horn in the
STRIP focal plane.

You should initialize `Horn` objects via the `InstrumentDB` constructor, which
loads their definition from a STRIP instrument database in YAML format.

Field              | Type           | Meaning
:----------------- |:-------------- |:------------------------------------------------------------
`name`             | String         | Name of the horn, e.g., `I0`
`id`               | Int            | Unique number of the horn, starting from 1
`polid`            | Int            | Unique ID of the polarimeter associated with the horn
`polarizerid`      | Int            | Unique ID of the polarizer+OMT associated with the horn
`moduleid`         | Int            | Number of the horn within the module, from 0 to 6
`color`            | String         | Name of the color associated with the module
`orientation`      | Array{Float64} | 3D vector containing the orientation of the horn in the sky
`fwhm_x_deg`       | Float64        | FWHM of the beam along the X axis, in degrees
`fwhm_y_deg`       | Float64        | FWHM of the beam along the Y axis, in degrees
`main_spillover`   | Float64        | Main reflector spillover
`sub_spillover`    | Float64        | Sub-reflector spillover
`xpd_db`           | Float64        | Cross-polarization, in dB
`directivity_dbi`  | Float64        | Directivity, in dBi
`ellipticity`      | Float64        | Ellipticity
"""
struct Horn
    name::String
    id::Int
    polid::Int
    polarizerid::Int
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

function Base.show(io::IO, horn::Horn)
    if get(io, :compact, false)
        print(io, "Horn $(horn.name), module $(horn.moduleid) ($(horn.color))")
    else
        @printf(io, """
            Horn %s, module %d (%s):
                Orientation: [%.4f, %.4f, %.4f]
                FWHM (X/Y): %.4f × %.4f °
                Spillover: %f (main), %f (sub)
                Cross-polarization: %.2f dB
                Directivity: %.2f dBi
                Ellipticity: %.4f""",
            horn.name, horn.moduleid, horn.color,
            horn.orientation[1], horn.orientation[2], horn.orientation[3],
            horn.fwhm_x_deg, horn.fwhm_y_deg, 
            horn.main_spillover, horn.sub_spillover,
            horn.xpd_db,
            horn.directivity_dbi,
            horn.ellipticity)
    end
end

@doc raw"""
    BandshapeInfo

Information about the spectral band response of a polarimeter.

Field                      | Type             | Meaning
:------------------------- |:---------------- |:-------------------------------------------------------
`center_frequency_hz`      | Float64          | Estimate for the center frequency, in Hz
`center_frequency_err_hz`  | Float64          | Estimated error on the center frequency, in Hz
`bandwidth_hz`             | Float64          | Estimated bandwidth, in Hz
`bandwidth_err_hz`         | Float64          | Estimated error on the bandwidth, in Hz
`lowest_frequency_hz`      | Float64          | Lowest frequency of the bandshape in `response`, in Hz
`highest_frequency_hz`     | Float64          | Highest frequency of the bandshape in `response`, in Hz
`num_of_frequencies`       | Int              | Number of samples in `response`
`bandshape`                | Array{Float64,1} | Profile of the bandshape (pure numbers)
`bandshape_error`          | Array{Float64,1} | Estimated error on the profile of the bandshape 
`test_id`                  | Array{Int,1}     | ID of the unit-level test used to characterize the bandshape
`analysis_id`              | Int              | ID of the unit-level analysis used to characterize the bandshape

You can plot a `BandshapeInfo` object by importing `Plots` and using `plot`:

```julia
db = InstrumentDB()
plot(bandpass(db, "I0"), show_error = true)
```

The following keywords are recognized in the call to `plot`:

- `show_error` (default: `true`): include an error bar.
- `show_centerfreq` (default: `false`): include a vertical bar showing
  the position of the center frequency

"""
struct BandshapeInfo
    center_frequency_hz::Float64
    center_frequency_err_hz::Float64
    bandwidth_hz::Float64
    bandwidth_err_hz::Float64
    lowest_frequency_hz::Float64
    highest_frequency_hz::Float64
    num_of_frequencies::Int
    bandshape::Array{Float64,1}
    bandshape_error::Array{Float64,1}
    test_id::Array{Int,1}
    analysis_id::Int
end

function Base.show(io::IO, band::BandshapeInfo)
    if get(io, :compact, false)
        @printf(io, "BandshapeInfo(ν0=%.2f GHz, Δν=%.2f GHz)",
                band.center_frequency_hz * 1e-9,
                band.bandwidth_hz * 1e-9)
    else
        @printf(io, """
            Bandshape:
                Center frequency: %.2f ± %.2f GHz
                Bandwidth: %.2f ± %.2f GHz
                Frequency range: [%.2f, %.2f] GHz (%d points)
                Test ID: [%s]
                Analysis ID: %d""",
            band.center_frequency_hz * 1e-9, band.center_frequency_err_hz * 1e-9,
            band.bandwidth_hz * 1e-9, band.bandwidth_err_hz * 1e-9,
            band.lowest_frequency_hz * 1e-9, band.highest_frequency_hz * 1e-9, band.num_of_frequencies,
            join(["$x" for x in band.test_id], ", "),
            band.analysis_id)
    end
end

RecipesBase.@recipe function plot(band::BandshapeInfo; show_error = true, show_centerfreq = false)
    seriestype --> :path
    xguide --> "Frequency [GHz]"

    if show_centerfreq
        RecipesBase.@series begin
            seriestype --> :path
            color --> :gray
            label --> ""
            let centerfreq = band.center_frequency_hz * 1e-9
                [centerfreq, centerfreq], [0, 1]
            end
        end
    end
    
    if show_error
        ribbon --> (band.bandshape_error, band.bandshape_error)
        fillalpha --> 0.3
    end
    
    ν, profile = bandshape(band)
    ν .* 1e-9, profile
end


@doc raw"""
    BandshapeInfo()

Initialize a BandshapeInfo object with all values set to zero.
"""
BandshapeInfo() = BandshapeInfo(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, Float64[], 0, 0)

function bandshape(bandinfo::BandshapeInfo)
    ν = range(
        bandinfo.lowest_frequency_hz,
        stop = bandinfo.highest_frequency_hz,
        length = bandinfo.num_of_frequencies,
    )

    @assert length(ν) == length(bandinfo.bandshape)
    (ν, bandinfo.bandshape, bandinfo.bandshape_error)
end

@doc raw"""
    bandshape(bandinfo::BandshapeInfo) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}
    bandshape(db::InstrumentDB, polid::Integer) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}
    bandshape(db::InstrumentDB, horn_name::AbstractString) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}

Return a pair `(ν_hz, B, Berr)` containing the shape of the bandpass
in `bandinfo` (first form), or the bandpass taken from the instrument
database (second and third form). The two elements of the tuple
`(ν_hz, B)` are two arrays of the same length containing the
frequencies (in Hz) and the bandpass response at the same frequency
(pure number), and they are suitable to be plotted, like in the
following example:

```julia
db = InstrumentDB()
x, y, err = bandshape(db, "G2")
plot(x, y, ribbon=(err, err))   # Plot the bandpass and the error bar
```

However, it is easier just to use `plot` on a [`BandshapeInfo`](@ref)
object.

"""
bandshape

@doc raw"""
    SpectrumInfo

Information about the noise spectrum of the output of a polarimeter.

Field                | Type     | Meaning
:------------------- |:-------- |:-------------------------------------------------------
`slope_i`            | Float64  | The slope ($\alpha$) of the 1/f component of the noise in the I signal
`slope_i_err`        | Float64  | Error associated with the value of `slope_i`
`slope_q`            | Float64  | Same as `slope_i`, but for the Q signal
`slope_q_err`        | Float64  | Error associated with the value of `slope_q`
`slope_u`            | Float64  | Same as `slope_i`, but for the U signal
`slope_u_err`        | Float64  | Error associated with the value of `slope_u`
`fknee_i_hz`         | Float64  | Knee frequency of the I signal, in Hz
`fknee_i_err_hz`     | Float64  | Error associated with the value of `fknee_i_hz`
`fknee_q_hz`         | Float64  | Knee frequency of the Q signal, in Hz
`fknee_q_err_hz`     | Float64  | Error associated with the value of `fknee_q_hz`
`fknee_u_hz`         | Float64  | Knee frequency of the U signal, in Hz
`fknee_u_err_hz`     | Float64  | Error associated with the value of `fknee_u_hz`
`wn_i_k2_hz`         | Float64  | White noise level for the I signal, in K^2 Hz
`wn_i_err_k2_hz`     | Float64  | Error associated with the value of `wn_i_k2_hz`
`wn_q_k2_hz`         | Float64  | White noise level for the Q signal, in K^2 Hz
`wn_q_err_k2_hz`     | Float64  | Error associated with the value of `wn_q_k2_hz`
`wn_u_k2_hz`         | Float64  | White noise level for the U signal, in K^2 Hz
`wn_u_err_k2_hz`     | Float64  | Error associated with the value of `wn_u_k2_hz`
`load_temperature_k` | Float64  | System brightness temperature used during the tests (in K)
`test_id`            | Int      | ID of the unit-level test used to characterize the bandshape
`analysis_id`        | Int      | ID of the unit-level analysis used to characterize the bandshape

You can quickly plot the theoretical shape of the noise power spectrum
using `plot` on a `SpectrumInfo` object.

"""
struct SpectrumInfo
    slope_i::Float64
    slope_q::Float64
    slope_u::Float64
    slope_i_err::Float64
    slope_q_err::Float64
    slope_u_err::Float64
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
    load_temperature_k::Float64
    test_id::Int
    analysis_id::Int
end

function Base.show(io::IO, spec::SpectrumInfo)
    if get(io, :compact, false)
        @printf(io, "SpectrumInfo(α=[%.4f, %.4f, %.4f], fknee=[%.3f, %.3f, %.3f], wn=[%.3f, %.3f, %.3f])",
                spec.slope_i, spec.slope_q, spec.slope_u,
                spec.fknee_i_hz, spec.fknee_q_hz, spec.fknee_u_hz,
                spec.wn_i_k2_hz, spec.wn_q_k2_hz, spec.wn_u_k2_hz)
    else
        @printf(io, """
            Noise spectrum:
                Slope: I = %.4f ± %.4f, Q = %.4f ± %.4f, U = %.4f ± %.4f
                Knee frequency: I = %.1f ± %.1f mHz, Q = %.1f ± %.1f mHz, U = %.1f ± %.1f mHz
                White noise: Q = %.1f ± %.1f mK^2 Hz, U = %.1f ± %.1f mK^2 Hz
                System brightness temperature: %.1f K
                Test ID: %d
                Analysis ID: %d""",
            spec.slope_i, spec.slope_i_err,
            spec.slope_q, spec.slope_q_err,
            spec.slope_u, spec.slope_u_err,
            spec.fknee_i_hz * 1e3, spec.fknee_i_err_hz * 1e3,
            spec.fknee_q_hz * 1e3, spec.fknee_q_err_hz * 1e3,
            spec.fknee_u_hz * 1e3, spec.fknee_u_err_hz * 1e3,
            spec.wn_q_k2_hz * 1e6, spec.wn_q_err_k2_hz * 1e6,
            spec.wn_u_k2_hz * 1e6, spec.wn_u_err_k2_hz * 1e6,
            spec.load_temperature_k,
            spec.test_id,
            spec.analysis_id)
    end
end

RecipesBase.@recipe function plot(spec::SpectrumInfo)
    model = (freq, α, fknee, wn) -> wn * (1 + fknee/freq)^α
    fknee = Float64[]
    spec.fknee_q_hz > 0 && push!(fknee, spec.fknee_q_hz)
    spec.fknee_u_hz > 0 && push!(fknee, spec.fknee_u_hz)

    min_freq, max_freq = if length(fknee) ≥ 1
        1e-4 * minimum(fknee), 1e+4 * maximum(fknee)
    else
        0.01, 1.0
    end
    freqs = 10 .^ range(log10(min_freq), log10(max_freq), length = 30)

    scale --> :log10
    
    RecipesBase.@series begin
        label --> "Q"
        (freqs, model.(freqs, Ref(spec.slope_q), Ref(spec.fknee_q_hz), Ref(spec.wn_q_k2_hz)))
    end

    label --> "U"
    xlabel --> "Frequency [Hz]"
    ylabel --> "Power [K^2/Hz]"
    (freqs, model.(freqs, Ref(spec.slope_u), Ref(spec.fknee_u_hz), Ref(spec.wn_u_k2_hz)))
end

@doc raw"""
    SpectrumInfo()

Initialize a SpectrumInfo object with all values set to zero.
"""
SpectrumInfo() = SpectrumInfo(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)

@doc raw"""
    NoiseTemperatureInfo

Information about the noise temperature of a polarimeter. This structure is used
for the field `tnoise` of the `Detector` struct.
    
Field          | Type             | Meaning
:------------- |:---------------- |:------------------------------------------------------------------------------
`tnoise_k`     | Float64          | Noise temperature computed from `tnoise_values_k`, in K
`tnoise_err_k` | Float64          | Error associated with `tnoise_k`, computed from `tnoise_values_k`
`test_ids`     | Array{Int,1}     | List of unit-level test IDs used to estimate the noise temperature
`analysis_ids` | Array{Int,1}     | List of unit-level analysis report IDs used to estimate the noise temperature
`values_k`     | Array{Float64,1} | List of noise temperatures estimated from the tests
"""
struct NoiseTemperatureInfo
    tnoise_k::Float64
    tnoise_err_k::Float64
    test_ids::Array{Int,1}
    analysis_ids::Array{Int,1}
    values_k::Array{Float64,1}
end

function Base.show(io::IO, tnoise::NoiseTemperatureInfo)
    if get(io, :compact, false)
        @printf(io, "NoiseTemperatureInfo(tnoise=%.1f)",
                tnoise.tnoise_k)
    else
        @printf(io, """
            Noise temperature:
                Tnoise: %.1f ± %.1f K
                Estimates: [%s] K
                Test IDs: [%s]
                Analysis IDs: [%s]""",
            tnoise.tnoise_k, tnoise.tnoise_err_k,
            join([@sprintf("%.1f", x) for x in tnoise.values_k], ", "),
            join(["$x" for x in tnoise.test_ids], ", "),
            join(["$x" for x in tnoise.analysis_ids], ", ")
            )
    end
end

@doc raw"""
    NoiseTemperatureInfo()

Initialize a NoiseTemperatureInfo object with all values set to zero.
"""
NoiseTemperatureInfo() = NoiseTemperatureInfo(0.0, 0.0, Int[], Int[], Float64[])

@doc raw"""
Information about a STRIP detector

This structure holds information about a STRIP polarimeter.

You should initialize `Detector` objects via the `InstrumentDB` constructor,
which loads their definition from a local STRIP instrument database.

Field       | Type                 | Meaning
:---------- |:-------------------- |:-----------------------------------------------------------------
`id`        | Int                  | Integer ID of the polarimeter, e.g., `2` for `STRIP02`
`name`      | String               | Full name of the polarimeter, e.g., `STRIP02`
`band`      | String               | Band: it can either be `Q` or `W`
`bandshape` | BandshapeInfo        | Information about the bandpass response
`spectrum`  | SpectrumInfo         | Information about the noise spectrum (white noise and 1/f noise)
`tnoise`    | NoiseTemperatureInfo | Information about the noise temperature
"""
struct Detector
    id::Int
    name::String
    band::String
    
    bandshape::BandshapeInfo
    spectrum::SpectrumInfo
    tnoise::NoiseTemperatureInfo
end

function Base.show(io::IO, det::Detector)
    if get(io, :compat, false)
        print(io, "Detector($(det.name), $(det.band) band)")
    else
        @printf(io, """
            Detector %s (%s band):
                Center frequency: %.2f ± %.2f GHz
                Bandwidth: %.2f ± %.2f GHz
                Noise temperature: %.1f ± %.2f K
                Knee frequency: %.1f ± %.1f mHz (Q), %.1f ± %.1f mHz (U)""",
            det.name, det.band,
            det.bandshape.center_frequency_hz * 1e-9, det.bandshape.center_frequency_err_hz * 1e-9,
            det.bandshape.bandwidth_hz * 1e-9, det.bandshape.bandwidth_err_hz * 1e-9,
            det.tnoise.tnoise_k, det.tnoise.tnoise_err_k,
            det.spectrum.fknee_q_hz * 1e3, det.spectrum.fknee_q_err_hz * 1e3,
            det.spectrum.fknee_u_hz * 1e3, det.spectrum.fknee_u_err_hz * 1e3)
    end
end

@doc raw"""
STRIP instrument database

The "database" contains information about feed horns and polarimeters:

- The field `focalplane` is a dictionary (mapping) associating the string
  identifying a horn (e.g., `I0`) with a [`Horn`](@ref) structure;

- The field `detectors` is a dictionary associating the ID of the
  polarimeter (e.g., 2 stands for `STRIP02`) with a [`Detector`](@ref) structure.

You should usually create an object of this kind using the default constructor,
which parses a set of YAML files containing the real parameters of the
instrument.

# Examples

```jldoctest
julia> db = InstrumentDB();

julia> print("Number of horns in the database: $(length(keys(db.focalplane)))")
Number of horns in the database: 55

julia> print("Number of polarimeters in the database: $(length(keys(db.detectors)))")
Number of polarimeters in the database: 66
```

# Visualization

You can produce a table describing the contents of the instrument
database using `show` and passing `text/markdown` as MIME type:

````julia
db = InstrumentDB()
show(stdout, MIME("text/markdown"), db)
````

The table can be converted to other formats (HTML, LaTeX, Microsoft
Word, …) using commonly-available tools, e.g.,
[Pandoc](https://pandoc.org/).

"""
struct InstrumentDB
    focalplane::Dict{String,Horn}
    detectors::Dict{Int,Detector}
end

function Base.show(io::IO, db::InstrumentDB)
    @printf(io, "InstrumentDB(%d horns, %d detectors)",
            length(keys(db.focalplane)),
            length(keys(db.detectors)))
end

function Base.show(io::IO, ::MIME"text/markdown", db::InstrumentDB)
    print(io, """# Horns

| Horn | Band | Polarimeter | FWHM (x) [deg] | FWHM (y) [deg] |
|------|------|-------------|----------------|----------------|
""")

    # Present the horns in a well-defined order (alphabetic order
    # would put W-band before the Y module).
    for mod in ["I", "V", "B", "G", "Y", "O", "R", "W"]
        for horn in 0:6
            cur_key = "$mod$horn"
            cur_key in keys(db.focalplane) || continue
            cur_horn = db.focalplane[cur_key]
            @printf(
                "| %s | %s | %s | %.2f | %.2f |\n",
                cur_key,
                db.detectors[cur_horn.polid].name,
                db.detectors[cur_horn.polid].band,
                cur_horn.fwhm_x_deg,
                cur_horn.fwhm_y_deg,
            )
        end
    end

    print(io, """

# Polarimeters

| Name | Band | Noise temperature [K] | Center Frequency [GHz] | Bandpass [GHz] |
|------|------|-----------------------|------------------------|----------------|
""")

    det_keys = sort(keys(db.detectors) |> collect)
    for cur_key in det_keys
        cur_horn = db.detectors[cur_key]
        @printf(
            "| %s | %s | %.1f | %.2f | %.2f |\n",
            cur_horn.name,
            cur_horn.band,
            cur_horn.tnoise.tnoise_k,
            cur_horn.bandshape.center_frequency_hz * 1e-9,
            cur_horn.bandshape.bandwidth_hz * 1e-9,
        )
    end
        
end

@doc raw"""
    parsefpdict(fpdict)

Return a dictionary associating an horn name (e.g., `I0`) to a `Horn`
object containing information about some horn in the STRIP focal
plane. The information are parsed from `fpdict`, which should be a
dictionary loaded from a YAML file. The default YAML file to be used
is located in the folder returned by [`defaultdbfolder`](@ref) and is
usually named `strip_focal_plane.yaml`

"""
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
                               value["polarizer_id"],
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
        get(banddict, "bandshape", Float64[]),
        get(banddict, "bandshape_error", Float64[]),
        get(banddict, "test_id", Int[]),
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
        get(specdict, "load_average_temperature_k", 20.0),
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

@doc raw"""
    parsedetdict(detdict)

Return a dictionary associating an integer number to a `Detector`
object containing information about the STRIP detector with the
corresponding number. The information are parsed from `detdict`, which
should be a dictionary loaded from a YAML file. The default YAML file
to be used is located in the folder returned by
[`defaultdbfolder`](@ref) and is usually named `strip_detectors.yaml`

"""
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

@doc raw"""
    InstrumentDB(dbpath::AbstractString)

Load the STRIP instrument database from the specified path. Return an
instance of a InstrumentDB object.

"""
function InstrumentDB(dbpath::AbstractString)
    focalplanedict = open(joinpath(dbpath, "strip_focal_plane.yaml")) do f
        YAML.load(f)
    end
    
    detectordict = open(joinpath(dbpath, "strip_detectors.yaml")) do f
        YAML.load(f)
    end
    
    InstrumentDB(parsefpdict(focalplanedict), parsedetdict(detectordict))
end

@doc raw"""
    defaultdbfolder()

Return a string containing the (local) full path to the YAML files
containing the reference instrument DB.

"""
defaultdbfolder() = joinpath(dirname(pathof(Stripeline)), "..", "instrumentdb")

@doc raw"""
    InstrumentDB()

Load the STRIP instrument database from the directory returned by
[`defaultdbfolder`](@ref). Return an instance of a InstrumentDB object.
"""
InstrumentDB() = InstrumentDB(defaultdbfolder())

################################################################################

detector(db::InstrumentDB, polid::Integer) = db.detectors[polid]
detector(db::InstrumentDB, horn_name::AbstractString) = detector(db, db.focalplane[horn_name].polid)

@doc raw"""
    detector(db::InstrumentDB, polid::Integer) -> Detector
    detector(db::InstrumentDB, horn_name::AbstractString) -> Detector

Return a [`Detector`](@ref) structure, taken from the instrument
database. If the form with `polid` is used, `polid` is the progressive
number of the polarimeter; e.g., for `STRIP02`, `polid == 2`. In the
second form, you pass the string identifying the horn on the focal
plane, e.g., `I0`, `W3`, etc.

```julia
db = InstrumentDB()
pol1 = detector(db, 16)   # Get information about STRIP16
pol2 = detector(db, "V4") # Get information about the detector connected to horn V4
```

"""
detector

# These functions are already documented above
bandshape(db::InstrumentDB, polid::Integer) = bandshape(bandpass(db, polid))
bandshape(db::InstrumentDB, horn_name::AbstractString) = bandshape(bandpass(db, horn_name))

bandpass(db::InstrumentDB, polid::Integer) = detector(db, polid).bandshape
bandpass(db::InstrumentDB, horn_name::AbstractString) = bandpass(db, db.focalplane[horn_name].polid)

@doc raw"""
    bandpass(db::InstrumentDB, polid::Integer) -> BandshapeInfo
    bandpass(db::InstrumentDB, horn_name::AbstractString) -> BandshapeInfo

Return a pair `(ν_hz, B)` containing the bandpass `B` for the horn
with the specified ID (`polid`) or associated to some horn
(`horn_name`). To understand how `polid` and `horn_name` work, see the
documentation for [`detector`](@ref).

The two elements of the tuple `(ν_hz, B)` are two arrays of the same
length containing the frequencies (in Hz) and the bandpass response at
the same frequency (pure number).

```julia
db = InstrumentDB()
x, y = bandpass(db, "G2")
plot(x, y)   # Plot the bandpass
```

See also [`bandshape`](@ref).

"""
bandpass

spectrum(db::InstrumentDB, polid::Integer) = detector(db, polid).spectrum
spectrum(db::InstrumentDB, horn_name::AbstractString) = spectrum(db, db.focalplane[horn_name].polid)

@doc raw"""
    spectrum(db::InstrumentDB, polid::Integer) -> SpectrumInfo
    spectrum(db::InstrumentDB, horn_name::AbstractString) -> SpectrumInfo

Return a [`SpectrumInfo`](@ref) object, taken from the instrument
database. The meaning of the parameters `polid` and `horn_name` is
explained in the documentation for [`detector`](@ref).

"""
spectrum

tnoise(db::InstrumentDB, polid::Integer) = detector(db, polid).tnoise
tnoise(db::InstrumentDB, horn_name::AbstractString) = tnoise(db, db.focalplane[horn_name].polid)

@doc raw"""
    tnoise(db::InstrumentDB, polid::Integer) -> NoiseTemperatureInfo
    tnoise(db::InstrumentDB, horn_name::AbstractString) -> NoiseTemperatureInfo

Return a [`NoiseTemperatureInfo`](@ref) object, taken from the
instrument database. The meaning of the parameters `polid` and
`horn_name` is explained in the documentation for [`detector`](@ref).

"""
tnoise

function fknee_hz(db::InstrumentDB, polid::Integer; tsys_k = missing)
    specinfo = spectrum(db, polid)

    tsys_k === missing && return (specinfo.fknee_i_hz, specinfo.fknee_q_hz, specinfo.fknee_u_hz)

    # Correct the value of fknee depending on the ratio between the system
    # temperature now and the temperature used during the characterization of the
    # polarimeter

    (α_i, α_q, α_u) = (specinfo.slope_i, specinfo.slope_q, specinfo.slope_u)
    load_ratio = specinfo.load_temperature_k / tsys_k

    (specinfo.fknee_i_hz * load_ratio^(1 / α_i),
        specinfo.fknee_q_hz * load_ratio^(1 / α_q),
        specinfo.fknee_u_hz * load_ratio^(1 / α_u),)
end

fknee_hz(
    db::InstrumentDB,
    horn_name::AbstractString;
    tsys_k = missing,
) = fknee_hz(db, db.focalplane[horn_name].polid; tsys_k = tsys_k)

@doc raw"""
    fknee_hz(db::InstrumentDB, polid::Integer; tsys_k = missing) -> Tuple{Float64, Float64, Float64}
    fknee_hz(db::InstrumentDB, horn_name::AbstractString; tsys_k = missing) -> Tuple{Float64, Float64, Float64}

Return the knee frequency for the selected detector, taken from the
instrument database. The meaning of the parameters `polid` and
`horn_name` is explained in the documentation for [`detector`](@ref).

If `tsys_k` is specified, the system temperature is rescaled to the
desired temperature of the load feeding the polarimeter, so that the
1/f component of the noise remains unchanged but the white noise
plateau raises/lowers by an appropriate amount. Otherwise, the
function returns the raw frequency taken from the instrument database.

"""
fknee_hz

################################################################################

@doc raw"""
    sensitivity_tant(db::InstrumentDB, load_tant; modules = Set([0, 1, 2, 3, 4, 5, 6]))

Calculate the white-noise sensitivity of an array of detectors,
measured in K⋅√s, given some antenna temperature for the load. The
result takes in account only those horns belonging to the modules
listed in the keyword `modules` (the W-band horns belong to module
`-1`). By default, only the Q-band modules are considered.

The result assumes the radiometer equation: ``σ⋅√τ =
\frac{T_{sys}}{√2β}``, where `T_{sys}` is the system temperature, `β`
is the bandwidth, and `τ` is the acquisition time. The factor 2 comes
from the way Strip polarimeters operate. The system temperature is
assumed to be the noise temperature of each detector, plus the term
`load_tant`, which should take into account all the other sources of
power entering the system (e.g., telescope, atmosphere, etc.). The
term `load_tant` should be expressed as an antenna temperature.

"""
function sensitivity_tant(db::InstrumentDB, load_tant; modules = Set([0, 1, 2, 3, 4, 5, 6]))
    result = 0.0
    num = 0

    for (name, horn) in db.focalplane
        (horn.moduleid ∈ modules) || continue

        polarimeter = db.detectors[horn.polid]
        tsys = load_tant + polarimeter.tnoise.tnoise_k
        result += 2 * polarimeter.bandshape.bandwidth_hz / (tsys^2)
        num += 1
    end

    (1.0 / sqrt(result), num)
end

const kb = 1.3806503e-23  # Boltzmann constant
const hplanck = 6.62607015e-34 # Planck constant

@doc raw"""
    t_to_trj(temperature_k, nu_hz)

Convert a thermodynamic temperature (in K) into a Rayleigh-Jeans temperature, given some
specified frequency `nu_hz` (in Hz).

See also [`trj_to_t`](@ref) for the inverse transformation.
"""
t_to_trj(temperature_k, nu_hz) = hplanck * nu_hz / kb / (exp(hplanck * nu_hz / (kb * temperature_k)) - 1)

function bisect(fn, range)
    # Plain old bisection method, with a fixed precision

    min_x, max_x = range
    threshold = 1e-9

    while true
        mid_x = (min_x + max_x) / 2
        abs(max_x - min_x) < threshold && return mid_x

        min_y, mid_y, max_y = (fn(x) for x in (min_x, mid_x, max_x))
        if sign(min_y * mid_y) < 0
            max_x = mid_x
        else
            min_x = mid_x
        end
    end
end


@doc raw"""
    trj_to_t(temperature_k, nu_hz)

Convert a Rayleigh-Jeans temperature (in K) into a thermodynamic temperature, given some
specified frequency `nu_hz` (in Hz).

See also [`t_to_trj`](@ref) for the inverse transformation.

"""
function trj_to_t(temperature_k, nu_hz)
    bisect(x->t_to_trj(x, nu_hz) - temperature_k,
        (temperature_k / 10.0, temperature_k * 10.0),
    )
end

@doc raw"""
    deltat_to_deltatrj(temperature_k, deltat_k, nu_hz)

Convert a small temperature fluctuation `deltat_k` around temperature
`temperature_k` from thermodynamic temperature to Rayleigh-Jeans (RJ)
temperature. This function can be used to convert sensitivities expressed as
thermodynamic temperatures in RJ sensitivities.

See also `deltatrj_to_deltat` for the inverse function.
"""
function deltat_to_deltatrj(temperature_k, deltat_k, nu_hz)
    # The formula we are using here is ∂T_RJ/∂T × δT, but we are
    # using the analytical formula for the derivative
    
    exponential = exp(hplanck * nu_hz / (kb * temperature_k))
    (hplanck * nu_hz / (temperature_k * kb))^2 * exponential / (exponential - 1)^2 * deltat_k
end

@doc raw"""
    deltat_to_deltatrj(temperature_k, deltat_k, nu_hz)

Convert a small temperature fluctuation `deltat_k` around temperature
`temperature_k` from Rayleigh-Jeans (RJ) temperature to thermodynamic
temperature. This function can be used to convert sensitivities expressed as
RJ temperatures in thermodynamic sensitivities.

See also `deltatrj_to_deltat` for the inverse function.
"""
function deltatrj_to_deltat(temperature_k, deltat_k, nu_hz)
    tcmb_k = trj_to_t(temperature_k, nu_hz)

    exponential = exp(hplanck * nu_hz / (kb * tcmb_k))
    (tcmb_k * kb / (hplanck * nu_hz))^2 * (exponential - 1)^2 / exponential * deltat_k
end
