var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "#Stripeline-User\'s-Manual-1",
    "page": "Introduction",
    "title": "Stripeline User\'s Manual",
    "category": "section",
    "text": "An implementation of a simulation/data analysis pipeline for the LSPE/STRIP instrument.To install Stripeline, start Julia and type the following command:using Pkg\nPkg.add(\"https://github.com/lspestrip/Stripeline\")To run the test suite, type the following command:using Pkg; Pkg.test(\"Stripeline\")In this manual, we will often assume that Stripeline has been imported using the following commands:import Stripeline\nconst Sl = StripelineIn this way, we can call functions like genpointings using the syntax Sl.genpointings, instead of the longer Stripeline.genpointings."
},

{
    "location": "#Documentation-1",
    "page": "Introduction",
    "title": "Documentation",
    "category": "section",
    "text": "The documentation was built using Documenter.jl.import Dates: now #hide\nprintln(\"Documentation built $(now()) with Julia $(VERSION).\") # hide"
},

{
    "location": "#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "basic/#",
    "page": "Basic functions",
    "title": "Basic functions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "basic/#Basic-functions-1",
    "page": "Basic functions",
    "title": "Basic functions",
    "category": "section",
    "text": ""
},

{
    "location": "basic/#Stripeline.sensitivity_tant",
    "page": "Basic functions",
    "title": "Stripeline.sensitivity_tant",
    "category": "function",
    "text": "sensitivity_tant(db::InstrumentDB, load_tant; modules = Set([0, 1, 2, 3, 4, 5, 6]))\n\nCalculate the white-noise sensitivity of an array of detectors, measured in K⋅√s, given some antenna temperature for the load. The result takes in account only those horns belonging to the modules listed in the keyword modules (the W-band horns belong to module -1). By default, only the Q-band modules are considered.\n\nThe result assumes the radiometer equation: στ = fracT_sys2β, where T_{sys} is the system temperature, β is the bandwidth, and τ is the acquisition time. The factor 2 comes from the way Strip polarimeters operate. The system temperature is assumed to be the noise temperature of each detector, plus the term load_tant, which should take into account all the other sources of power entering the system (e.g., telescope, atmosphere, etc.). The term load_tant should be expressed as an antenna temperature.\n\n\n\n\n\n"
},

{
    "location": "basic/#Stripeline.t_to_trj",
    "page": "Basic functions",
    "title": "Stripeline.t_to_trj",
    "category": "function",
    "text": "t_to_trj(temperature_k, nu_hz)\n\nConvert a thermodynamic temperature (in K) into a Rayleigh-Jeans temperature, given some specified frequency nu_hz (in Hz).\n\nSee also trj_to_t for the inverse transformation.\n\n\n\n\n\n"
},

{
    "location": "basic/#Stripeline.trj_to_t",
    "page": "Basic functions",
    "title": "Stripeline.trj_to_t",
    "category": "function",
    "text": "trj_to_t(temperature_k, nu_hz)\n\nConvert a Rayleigh-Jeans temperature (in K) into a thermodynamic temperature, given some specified frequency nu_hz (in Hz).\n\nSee also t_to_trj for the inverse transformation.\n\n\n\n\n\n"
},

{
    "location": "basic/#Stripeline.deltat_to_deltatrj",
    "page": "Basic functions",
    "title": "Stripeline.deltat_to_deltatrj",
    "category": "function",
    "text": "deltat_to_deltatrj(temperature_k, deltat_k, nu_hz)\n\nConvert a small temperature fluctuation deltat_k around temperature temperature_k from thermodynamic temperature to Rayleigh-Jeans (RJ) temperature. This function can be used to convert sensitivities expressed as thermodynamic temperatures in RJ sensitivities.\n\nSee also deltatrj_to_deltat for the inverse function.\n\n\n\n\n\n"
},

{
    "location": "basic/#Stripeline.deltatrj_to_deltat",
    "page": "Basic functions",
    "title": "Stripeline.deltatrj_to_deltat",
    "category": "function",
    "text": "deltat_to_deltatrj(temperature_k, deltat_k, nu_hz)\n\nConvert a small temperature fluctuation deltat_k around temperature temperature_k from Rayleigh-Jeans (RJ) temperature to thermodynamic temperature. This function can be used to convert sensitivities expressed as RJ temperatures in thermodynamic sensitivities.\n\nSee also deltatrj_to_deltat for the inverse function.\n\n\n\n\n\n"
},

{
    "location": "basic/#Measure-unit-conversion-1",
    "page": "Basic functions",
    "title": "Measure unit conversion",
    "category": "section",
    "text": "It is often useful to convert measurements between thermodynamic temperatures and Rayleigh-Jeans temperatures. The following functions implement this kind of conversion. Note that there are two families of functions:Functions that convert absolute measurements: t_to_trj, trj_to_t;\nFunctions that convert sensitivities (i.e., small fluctuations around an absolute value): deltat_to_deltatrj, deltatrj_to_deltatThe function sensitivity_tant computes the overall sensitivity of a set of polarimeters, using information from the Instrument database.sensitivity_tant\nt_to_trj\ntrj_to_t\ndeltat_to_deltatrj\ndeltatrj_to_deltat"
},

{
    "location": "instrumentdb/#",
    "page": "Instrument database",
    "title": "Instrument database",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "instrumentdb/#Instrument-database-1",
    "page": "Instrument database",
    "title": "Instrument database",
    "category": "section",
    "text": "InstrumentDB takes advantage of the structures Detector and Horn to retrieve information about feed horns and detectors from a YAML file. There are a set of YAML files containing the default configuration for the STRIP instrument in the repository."
},

{
    "location": "instrumentdb/#Quick-introduction-1",
    "page": "Instrument database",
    "title": "Quick introduction",
    "category": "section",
    "text": "The following example initializes an object of type InstrumentDB with the values referred to the standard STRIP instrument:using Stripeline; # hide\ndb = InstrumentDB()As db is a struct, its field can be accessed with the usual dot notation. The two fields in db are focalplane and detectors. They are both dictionaries, associating horn names to Horn objects and detectors IDs to Detector objects, respectively:db.focalplane[\"I0\"]\ndb.detectors[2]A number of high-level functions ease the access of the fields in a InstrumentDB object:detector returns a Detector structure, containing the details of a polarimeter;\nbandpass returns a BandshapeInfo structure, containing the shape of the bandpass of a detector;\nspectrum returns a SpectrumInfo\nfknee returns the knee frequency of the 1/f noise for the I, Q, and U signals, adapted to the brightness temperature of the load being observed by the detector;\ntnoise returns the noise temperature for the I, Q, and U components.The structure Detector uses three structures to organize its data in a hierarchical way:BandshapeInfo\nSpectrumInfo\nNoiseTemperatureInfoAll these structures know how to show themselves on the REPL:db.detectors[2].bandshape\ndb.detectors[2].spectrum\ndb.detectors[2].tnoiseFor more information about the fields in the structures listed above, as well as their meaning, keep reading."
},

{
    "location": "instrumentdb/#Stripeline.InstrumentDB",
    "page": "Instrument database",
    "title": "Stripeline.InstrumentDB",
    "category": "type",
    "text": "STRIP instrument database\n\nThe \"database\" contains information about feed horns and polarimeters:\n\nThe field focalplane is a dictionary (mapping) associating the string identifying a horn (e.g., I0) with a Horn structure;\nThe field detectors is a dictionary associating the ID of the polarimeter (e.g., 2 stands for STRIP02) with a Detector structure.\n\nYou should usually create an object of this kind using the default constructor, which parses a set of YAML files containing the real parameters of the instrument.\n\nExamples\n\njulia> db = InstrumentDB();\n\njulia> print(\"Number of horns in the database: $(length(keys(db.focalplane)))\")\nNumber of horns in the database: 55\n\njulia> print(\"Number of polarimeters in the database: $(length(keys(db.detectors)))\")\nNumber of polarimeters in the database: 66\n\nVisualization\n\nYou can produce a table describing the contents of the instrument database using show and passing text/markdown as MIME type:\n\ndb = InstrumentDB()\nshow(stdout, MIME(\"text/markdown\"), db)\n\nThe table can be converted to other formats (HTML, LaTeX, Microsoft Word, …) using commonly-available tools, e.g., Pandoc.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.Horn",
    "page": "Instrument database",
    "title": "Stripeline.Horn",
    "category": "type",
    "text": "Information about a STRIP horn\n\nThis structure holds a number of parameters relative to each feed horn in the STRIP focal plane.\n\nYou should initialize Horn objects via the InstrumentDB constructor, which loads their definition from a STRIP instrument database in YAML format.\n\nField Type Meaning\nname String Name of the horn, e.g., I0\nid Int Unique number of the horn, starting from 1\npolid Int Unique ID of the polarimeter associated with the horn\npolarizerid Int Unique ID of the polarizer+OMT associated with the horn\nmoduleid Int Number of the horn within the module, from 0 to 6\ncolor String Name of the color associated with the module\norientation Array{Float64} 3D vector containing the orientation of the horn in the sky\nfwhm_x_deg Float64 FWHM of the beam along the X axis, in degrees\nfwhm_y_deg Float64 FWHM of the beam along the Y axis, in degrees\nmain_spillover Float64 Main reflector spillover\nsub_spillover Float64 Sub-reflector spillover\nxpd_db Float64 Cross-polarization, in dB\ndirectivity_dbi Float64 Directivity, in dBi\nellipticity Float64 Ellipticity\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.Detector",
    "page": "Instrument database",
    "title": "Stripeline.Detector",
    "category": "type",
    "text": "Information about a STRIP detector\n\nThis structure holds information about a STRIP polarimeter.\n\nYou should initialize Detector objects via the InstrumentDB constructor, which loads their definition from a local STRIP instrument database.\n\nField Type Meaning\nid Int Integer ID of the polarimeter, e.g., 2 for STRIP02\nname String Full name of the polarimeter, e.g., STRIP02\nband String Band: it can either be Q or W\nbandshape BandshapeInfo Information about the bandpass response\nspectrum SpectrumInfo Information about the noise spectrum (white noise and 1/f noise)\ntnoise NoiseTemperatureInfo Information about the noise temperature\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.BandshapeInfo",
    "page": "Instrument database",
    "title": "Stripeline.BandshapeInfo",
    "category": "type",
    "text": "BandshapeInfo\n\nInformation about the spectral band response of a polarimeter.\n\nField Type Meaning\ncenter_frequency_hz Float64 Estimate for the center frequency, in Hz\ncenter_frequency_err_hz Float64 Estimated error on the center frequency, in Hz\nbandwidth_hz Float64 Estimated bandwidth, in Hz\nbandwidth_err_hz Float64 Estimated error on the bandwidth, in Hz\nlowest_frequency_hz Float64 Lowest frequency of the bandshape in response, in Hz\nhighest_frequency_hz Float64 Highest frequency of the bandshape in response, in Hz\nnum_of_frequencies Int Number of samples in response\nbandshape Array{Float64,1} Profile of the bandshape (pure numbers)\nbandshape_error Array{Float64,1} Estimated error on the profile of the bandshape\ntest_id Array{Int,1} ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\nYou can plot a BandshapeInfo object by importing Plots and using plot:\n\ndb = InstrumentDB()\nplot(bandpass(db, \"I0\"), show_error = true)\n\nThe following keywords are recognized in the call to plot:\n\nshow_error (default: true): include an error bar.\nshow_centerfreq (default: false): include a vertical bar showing the position of the center frequency\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.bandshape",
    "page": "Instrument database",
    "title": "Stripeline.bandshape",
    "category": "function",
    "text": "bandshape(bandinfo::BandshapeInfo) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}\nbandshape(db::InstrumentDB, polid::Integer) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}\nbandshape(db::InstrumentDB, horn_name::AbstractString) -> Tuple{Array{Float64, 1}, Array{Float64, 1}}\n\nReturn a pair (ν_hz, B, Berr) containing the shape of the bandpass in bandinfo (first form), or the bandpass taken from the instrument database (second and third form). The two elements of the tuple (ν_hz, B) are two arrays of the same length containing the frequencies (in Hz) and the bandpass response at the same frequency (pure number), and they are suitable to be plotted, like in the following example:\n\ndb = InstrumentDB()\nx, y, err = bandshape(db, \"G2\")\nplot(x, y, ribbon=(err, err))   # Plot the bandpass and the error bar\n\nHowever, it is easier just to use plot on a BandshapeInfo object.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.SpectrumInfo",
    "page": "Instrument database",
    "title": "Stripeline.SpectrumInfo",
    "category": "type",
    "text": "SpectrumInfo\n\nInformation about the noise spectrum of the output of a polarimeter.\n\nField Type Meaning\nslope_i Float64 The slope (alpha) of the 1/f component of the noise in the I signal\nslope_i_err Float64 Error associated with the value of slope_i\nslope_q Float64 Same as slope_i, but for the Q signal\nslope_q_err Float64 Error associated with the value of slope_q\nslope_u Float64 Same as slope_i, but for the U signal\nslope_u_err Float64 Error associated with the value of slope_u\nfknee_i_hz Float64 Knee frequency of the I signal, in Hz\nfknee_i_err_hz Float64 Error associated with the value of fknee_i_hz\nfknee_q_hz Float64 Knee frequency of the Q signal, in Hz\nfknee_q_err_hz Float64 Error associated with the value of fknee_q_hz\nfknee_u_hz Float64 Knee frequency of the U signal, in Hz\nfknee_u_err_hz Float64 Error associated with the value of fknee_u_hz\nwn_i_k2_hz Float64 White noise level for the I signal, in K^2 Hz\nwn_i_err_k2_hz Float64 Error associated with the value of wn_i_k2_hz\nwn_q_k2_hz Float64 White noise level for the Q signal, in K^2 Hz\nwn_q_err_k2_hz Float64 Error associated with the value of wn_q_k2_hz\nwn_u_k2_hz Float64 White noise level for the U signal, in K^2 Hz\nwn_u_err_k2_hz Float64 Error associated with the value of wn_u_k2_hz\nload_temperature_k Float64 System brightness temperature used during the tests (in K)\ntest_id Int ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\nYou can quickly plot the theoretical shape of the noise power spectrum using plot on a SpectrumInfo object.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.NoiseTemperatureInfo",
    "page": "Instrument database",
    "title": "Stripeline.NoiseTemperatureInfo",
    "category": "type",
    "text": "NoiseTemperatureInfo\n\nInformation about the noise temperature of a polarimeter. This structure is used for the field tnoise of the Detector struct.\n\nField Type Meaning\ntnoise_k Float64 Noise temperature computed from tnoise_values_k, in K\ntnoise_err_k Float64 Error associated with tnoise_k, computed from tnoise_values_k\ntest_ids Array{Int,1} List of unit-level test IDs used to estimate the noise temperature\nanalysis_ids Array{Int,1} List of unit-level analysis report IDs used to estimate the noise temperature\nvalues_k Array{Float64,1} List of noise temperatures estimated from the tests\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Structures-1",
    "page": "Instrument database",
    "title": "Structures",
    "category": "section",
    "text": "InstrumentDB\nHorn\nDetector\nBandshapeInfo\nbandshape\nSpectrumInfo\nNoiseTemperatureInfo"
},

{
    "location": "instrumentdb/#Stripeline.detector",
    "page": "Instrument database",
    "title": "Stripeline.detector",
    "category": "function",
    "text": "detector(db::InstrumentDB, polid::Integer) -> Detector\ndetector(db::InstrumentDB, horn_name::AbstractString) -> Detector\n\nReturn a Detector structure, taken from the instrument database. If the form with polid is used, polid is the progressive number of the polarimeter; e.g., for STRIP02, polid == 2. In the second form, you pass the string identifying the horn on the focal plane, e.g., I0, W3, etc.\n\ndb = InstrumentDB()\npol1 = detector(db, 16)   # Get information about STRIP16\npol2 = detector(db, \"V4\") # Get information about the detector connected to horn V4\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.bandpass",
    "page": "Instrument database",
    "title": "Stripeline.bandpass",
    "category": "function",
    "text": "bandpass(db::InstrumentDB, polid::Integer) -> BandshapeInfo\nbandpass(db::InstrumentDB, horn_name::AbstractString) -> BandshapeInfo\n\nReturn a pair (ν_hz, B) containing the bandpass B for the horn with the specified ID (polid) or associated to some horn (horn_name). To understand how polid and horn_name work, see the documentation for detector.\n\nThe two elements of the tuple (ν_hz, B) are two arrays of the same length containing the frequencies (in Hz) and the bandpass response at the same frequency (pure number).\n\ndb = InstrumentDB()\nx, y = bandpass(db, \"G2\")\nplot(x, y)   # Plot the bandpass\n\nSee also bandshape.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.spectrum",
    "page": "Instrument database",
    "title": "Stripeline.spectrum",
    "category": "function",
    "text": "spectrum(db::InstrumentDB, polid::Integer) -> SpectrumInfo\nspectrum(db::InstrumentDB, horn_name::AbstractString) -> SpectrumInfo\n\nReturn a SpectrumInfo object, taken from the instrument database. The meaning of the parameters polid and horn_name is explained in the documentation for detector.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.fknee_hz",
    "page": "Instrument database",
    "title": "Stripeline.fknee_hz",
    "category": "function",
    "text": "fknee_hz(db::InstrumentDB, polid::Integer; tsys_k = missing) -> Tuple{Float64, Float64, Float64}\nfknee_hz(db::InstrumentDB, horn_name::AbstractString; tsys_k = missing) -> Tuple{Float64, Float64, Float64}\n\nReturn the knee frequency for the selected detector, taken from the instrument database. The meaning of the parameters polid and horn_name is explained in the documentation for detector.\n\nIf tsys_k is specified, the system temperature is rescaled to the desired temperature of the load feeding the polarimeter, so that the 1/f component of the noise remains unchanged but the white noise plateau raises/lowers by an appropriate amount. Otherwise, the function returns the raw frequency taken from the instrument database.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.tnoise",
    "page": "Instrument database",
    "title": "Stripeline.tnoise",
    "category": "function",
    "text": "tnoise(db::InstrumentDB, polid::Integer) -> NoiseTemperatureInfo\ntnoise(db::InstrumentDB, horn_name::AbstractString) -> NoiseTemperatureInfo\n\nReturn a NoiseTemperatureInfo object, taken from the instrument database. The meaning of the parameters polid and horn_name is explained in the documentation for detector.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#High-level-access-functions-1",
    "page": "Instrument database",
    "title": "High-level access functions",
    "category": "section",
    "text": "detector\nbandpass\nspectrum\nfknee_hz\ntnoise"
},

{
    "location": "instrumentdb/#Stripeline.defaultdbfolder",
    "page": "Instrument database",
    "title": "Stripeline.defaultdbfolder",
    "category": "function",
    "text": "defaultdbfolder()\n\nReturn a string containing the (local) full path to the YAML files containing the reference instrument DB.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.parsefpdict",
    "page": "Instrument database",
    "title": "Stripeline.parsefpdict",
    "category": "function",
    "text": "parsefpdict(fpdict)\n\nReturn a dictionary associating an horn name (e.g., I0) to a Horn object containing information about some horn in the STRIP focal plane. The information are parsed from fpdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_focal_plane.yaml\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Stripeline.parsedetdict",
    "page": "Instrument database",
    "title": "Stripeline.parsedetdict",
    "category": "function",
    "text": "parsedetdict(detdict)\n\nReturn a dictionary associating an integer number to a Detector object containing information about the STRIP detector with the corresponding number. The information are parsed from detdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_detectors.yaml\n\n\n\n\n\n"
},

{
    "location": "instrumentdb/#Loading-custom-databases-1",
    "page": "Instrument database",
    "title": "Loading custom databases",
    "category": "section",
    "text": "It is not needed to load the default instrument database, as Stripeline provides a number of additional functions to build mock databases from dictionaries.defaultdbfolder\nparsefpdict\nparsedetdict"
},

{
    "location": "scanning/#",
    "page": "Scanning strategy",
    "title": "Scanning strategy",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend<iframe\n    src=\"https://www.google.com/maps/embed?pb=!1m14!1m8!1m3!1d28777729.760743093!2d-34.4392714819374!3d28.30115741125523!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x0%3A0x0!2zMjjCsDE4JzAwLjkiTiAxNsKwMzAnMzYuNCJX!5e0!3m2!1sen!2sit!4v1586614115629!5m2!1sen!2sit\" \n    width=\"600\"\n    height=\"450\"\n    frameborder=\"0\"\n    style=\"border:0;\"\n    allowfullscreen=\"\"\n    aria-hidden=\"false\"\n    tabindex=\"0\">\n</iframe>"
},

{
    "location": "scanning/#Stripeline.telescopetoground",
    "page": "Scanning strategy",
    "title": "Stripeline.telescopetoground",
    "category": "function",
    "text": "telescopetoground(wheelanglesfn, time_s)\n\nReturn a quaternion of type Quaternion{Float64} representing the coordinate transform from the focal plane to the ground of the telescope. The parameter wheelanglesfn must be a function which takes as input a time, time_s, in seconds, and it must return a 3-tuple containing the angles of the following motors, measured in radians:\n\nThe boresight motor (rotation around the z axis, counterclockwise)\nThe altitude motor (rotation around the y axis, counterclockwise)\nThe ground motor (rotation around the z axis, clockwise: N→E→S→W)\n\nExample\n\ntelescopetoground(3600.0) do\n    # Boresight motor keeps a constant angle equal to 0°\n    # Altitude motor remains at 20° from the Zenith\n    # Ground motor spins at 1 RPM\n    (0.0, deg2rad(20.0), timetorotang(time_s, 1))\nend\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.groundtoearth",
    "page": "Scanning strategy",
    "title": "Stripeline.groundtoearth",
    "category": "function",
    "text": "groundtoearth(groundq, time_s, latitude_deg; day_duration_s=86400.0)\n\nReturn a quaternion of type Quaternion{Float64} representing the coordinate transformation from the ground of the telescope to the Equatorial coordinate system. The parameter groundq must be a quaternion describing the coordinate transformation from the focal plane of the telescope to the ground. The parameter time_s must be a time in seconds, and latitude_deg is the latitude (in degrees, N is positive) of the location where the observation is made.\n\nThe keyword day_duration_s specifies the length of a sidereal day in seconds.\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.genpointings",
    "page": "Scanning strategy",
    "title": "Stripeline.genpointings",
    "category": "function",
    "text": "genpointings!(wheelanglesfn, beam_dir, timerange_s, dirs, psi; \n              polaxis = Float64[1.0, 0.0, 0.0],\n              latitude_deg = TENERIFE_LATITUDE_DEG, \n              ground = false)\ngenpointings(wheelanglesfn, beam_dir, timerange_s; \n             polaxis = Float64[1.0, 0.0, 0.0],\n             latitude_deg = TENERIFE_LATITUDE_DEG,\n             ground = false)\ngenpointings!(wheelanglesfn, beam_dir, timerange_s, t_start, dirs, psi;\n              polaxis = Float64[1.0, 0.0, 0.0],\n              latitude_deg = TENERIFE_LATITUDE_DEG, \n              longitude_deg = TENERIFE_LONGITUDE_DEG,\n              height_m = TENERIFE_HEIGHT_M,\n              precession = true,\n              nutation = true, \n              aberration = true, \n              refraction = true)\ngenpointings(wheelanglesfn, beam_dir, timerange_s, t_start;\n             polaxis=Float64[1.0, 0.0, 0.0],\n             latitude_deg=TENERIFE_LATITUDE_DEG, \n             longitude_deg=TENERIFE_LONGITUDE_DEG,\n             height_m=TENERIFE_HEIGHT_M,\n             precession = true,\n             nutation = true, \n             aberration = true, \n             refraction = true)\n\nGenerate a set of pointing directions for a STRIP detector. Each function is provided in two flavours: the ones ending with ! save the results in the last two parameters dirs and psi, while the others automatically allocate memory and return their results as a pair (dirs, psi).\n\nThe parameter wheelanglesfn must be a function which takes as input a time in seconds and returns a 3-tuple containing the angles (in radians) of the three motors:\n\nThe boresight motor (rotation around the z axis, counterclockwise)\nThe altitude motor (rotation around the y axis, counterclockwise)\nThe ground motor (rotation around the z axis, clockwise: N→E→S→W)\n\nThe meaning of the parameters/keywords is the following:\n\nbeam_dir specifies the pointing direction of the mean (the boresight is [0, 0, 1]). It must be normalized.\ntimerange_s is an enumerable type that specifies at which times (in seconds) pointings must be computed.\nt_start is a DateTime which tells the exact UTC date and time of the  observation.\nlatitude_deg is the latitude of the location where the observation is made (in degrees, North is positive).\nground is a Boolean: if true, the function will return a 4-tuple containing the colatitude and longitude measured in Equatorial coordinates (columns 1 and 2) and in ground coordinates (columns 3 and 4). If false, only the Equatorial coordinates are computed.\npolaxis is the polarization axis; it must be normalized.\nprecession: if true (the default), the Earth\'s precession is taken into account.\nnutation: if true (the default), the Earth\'s nutation is taken into account.\naberration: if true (the default), stellar aberration is taken into account.\nrefraction: if true, refraction corrections are taken into account. As these corrections are only valid for optical wavelengths, the default is false.\n\nReturn values\n\nIf t_start is not provided, the function genpointings returns a 2-tuple containing the sky directions (a N×2 array containing declination and right ascension, in Equatorial coordinates) and the polarization angle for each time step. The function genpointings! sets the values in the last two parameters dirs and psi.\n\nIf t_start is provided, the function genpointings returns a 2-tuple (4-tuple) containing the directions (a N×2 or Nx4 array containing the colatitude and the longitude) and the polarization angles at each time step; genpointings! works as above.\n\nExamples\n\nHere is an example using the form without t_start:\n\ndir, psi = genpointings([0, 0, 1], 0:0.1:1) do time_s\n    # Boresight motor keeps a constant angle equal to 0°\n    # Altitude motor remains at 20° from the Zenith\n    # Ground motor spins at 1 RPM\n    (0.0, deg2rad(20.0), timetorotang(time_s, 1))\nend\n\nAnd here is an example using t_start; unlike the previous example, we use a lambda function instead of the do...end notation.\n\nimport Dates\n\ndirs, psi = genpointings(time_s -> (0, deg2rad(20),\n                                    timetorotang(time_s, 1)),\n                         [0, 0, 1],\n                         0:0.1:1,\n                         Dates.DateTime(2022, 01, 01, 0, 0, 0),\n                         latitude_deg=10.0,\n                         longitude_deg=20.0,\n                         height_m = 1000) do time_s\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.timetorotang",
    "page": "Scanning strategy",
    "title": "Stripeline.timetorotang",
    "category": "function",
    "text": "timetorotang(time, rpm)\n\nConvert a time into a rotation angle, given the number of rotations per minute. The time should be expressed in seconds. The return value is in radians. time can either be a scalar or a vector.\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.northdir",
    "page": "Scanning strategy",
    "title": "Stripeline.northdir",
    "category": "function",
    "text": "northdir(θ, ϕ)\neastdir(θ, ϕ)\n\nCompute the North/East versor for a vector. The North for a vector v is -dv/dθ, as θ is the colatitude and moves along the meridian, and the East is dv/dϕ.\n\nExamples\n\njulia> northdir(π/2, 0) ≈ [0, 0, 1]\ntrue\njulia> eastdir(π/2, 0) ≈ [1, 0, 0]\ntrue\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.eastdir",
    "page": "Scanning strategy",
    "title": "Stripeline.eastdir",
    "category": "function",
    "text": "northdir(θ, ϕ)\neastdir(θ, ϕ)\n\nCompute the North/East versor for a vector. The North for a vector v is -dv/dθ, as θ is the colatitude and moves along the meridian, and the East is dv/dϕ.\n\nExamples\n\njulia> northdir(π/2, 0) ≈ [0, 0, 1]\ntrue\njulia> eastdir(π/2, 0) ≈ [1, 0, 0]\ntrue\n\n\n\n\n\n"
},

{
    "location": "scanning/#Stripeline.polarizationangle",
    "page": "Scanning strategy",
    "title": "Stripeline.polarizationangle",
    "category": "function",
    "text": "polarizationangle(northdir, eastdir, poldir)\n\nCalculate the polarization angle projected in the sky in IAU conventions. The parameters northdir and eastdir must be versors that point towards the North and East, respectively; poldir must be a versor that identify the polarization direction projected in the sky. The return value is in radians, and it is zero if the polarization angles points toward East, π/2 if it points toward North, etc.\n\nExamples\n\njulia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 1, 0])\n0.0\njulia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 0, 1]) |> rad2deg\n90.0\njulia> polarizationangle([0, 0, 1], [0, 1, 0], [0, 0, -1]) |> rad2deg\n-90.0\n\n\n\n\n\n"
},

{
    "location": "scanning/#Simulating-the-scanning-strategy-1",
    "page": "Scanning strategy",
    "title": "Simulating the scanning strategy",
    "category": "section",
    "text": "telescopetoground\ngroundtoearth\ngenpointings\ntimetorotang\nnorthdir\neastdir\npolarizationangle"
},

{
    "location": "acquisition/#",
    "page": "Data acquisition",
    "title": "Data acquisition",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "acquisition/#Simulating-data-acquisition-1",
    "page": "Data acquisition",
    "title": "Simulating data acquisition",
    "category": "section",
    "text": "An essential part of Strip polarimeters is the set of four Analogue-to-Digital Converters (ADC) that measure input voltages as 20-bit digital numbers. This process is called quantization, and it is a non-invertible operation that causes loss of information. Ideal ADCs perform a linear operation (modulo a rounding operation), but real-world components are never perfectly linear.Because of the fact that CMB experiments like Strip measure brightness temperatures, Stripeline models ADCs as devices that convert temperatures into ADUs, neglecting the fact that Strip polarimeters convert incoming fluxes into voltages.Stripeline offers a few functions to simulate the behaviour of an ADC. The simulation of ADC behaviour is useful to estimate the impact of quantization and non linearities. The following schema show how things work:(Image: )Function adc_response takes a temperature as input, and it produces the output that would be emitted by an ADC. The function adc_inv_response performs the reverse transformation: it converts a digital number back to a temperature. Function adc_filter combines the two functions: it takes a temperature as input, and it returns the temperature that has been measured by the ADC, including the effect of quantization and non linearities.An important difference between adc_response and adc_inv_response is the fact that adc_response considers non linearities, while adc_inv_response does not.ADC\noptimal_adc\nadc_response\nadc_inv_response\nadc_filter"
},

{
    "location": "mapmaking/#",
    "page": "Map making",
    "title": "Map making",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "mapmaking/#Stripeline.tod2map_mpi",
    "page": "Map making",
    "title": "Stripeline.tod2map_mpi",
    "category": "function",
    "text": "tod2map(pix_idx, tod, num_of_pixels; comm = nothing) -> binned_map\ntod2map(pix_idx, tod, num_of_pixels, data_properties; comm = nothing) -> binned_map\n\nThis function creates a binned map from the time-ordered data kept in the array tod, assuming that each sample is observing the pixel whose index is in pix_idx. The parameter num_of_pixels contains the number of pixels in the Healpix map, and it is used as an upper bound for the values in pix_idx. The parameter comm must be a MPI communicator, or nothing if you are not using MPI.\n\nIf comm is not nothing, the function is parallelized using MPI, and each process computes a map from its available data.  All partial maps are then combined together with MPI.allreduce. The function returns an array containing the binned map.\n\nIf Array of structures TodNoisePropertiesis passed to the function, the output binned map will be a weighted binned map. Each sample is weighted according to the inverse white noise variance sigma^2 of the corrispondingù polarimeter. In this way, the less noisy polarimeters will count more in the estimation of the map.\n\nRequirements\n\nThe length of the arrays pix_idx and tod must be the same\n\n\n\n\n\n"
},

{
    "location": "mapmaking/#Stripeline.baseline2map_mpi",
    "page": "Map making",
    "title": "Stripeline.baseline2map_mpi",
    "category": "function",
    "text": "baseline2map_mpi(pix_idx, baselines, baseline_lengths, num_of_pixels; comm = nothing)-> noise_map\nbaseline2map_mpi(pix_idx, baselines, num_of_pixels, data_properties; comm = nothing) -> noise_map\n\nThis is a MPI based function: each MPI process computes a map from its available data.  All partial maps are then combined together with MPI.allreduce. The function returns an array containing the binned map.\n\nIf Array of structures TodNoisePropertiesis passed to the function (instead of baseline_lengths) the output binned map will be a weighted binned map. Each sample will be weighted according to the inverse white noise variance sigma^2 of the corrisponding polarimeter. In this way, the less noisy polarimeters will count more in the estimation of the map.\n\nRequirements\n\nThe length of baselines and baseline_lengths must be the same;\nThe value sum(baseline_lengths) must be the same as the length of pix_idx.\n\n\n\n\n\n"
},

{
    "location": "mapmaking/#Stripeline.destripe",
    "page": "Map making",
    "title": "Stripeline.destripe",
    "category": "function",
    "text": "destripe(pix_idx, tod, num_of_pixels, data_properties, rank; comm = nothing, threshold = 1e-9, max_iter = 10000, save_baseline_history = false) -> DestripingResults\n\nThis MPI based function creates a map from a TOD and removes both 1/f and white noise, using the destriping technique.\n\nThe parameters passed to the function have the following meaning:\n\npix_idx: array containing the indices of the pixels visited by the instrument\ntod: the values measured by the polarimeters for each pixel (either I, Q, or U)\nnum_of_pixels: the number of pixels in the map to be  produced. This is used as an upper limit for the values in pix_idx\ndata_properties: an array of structures TodNoiseProperties,  holding information on each simulated polarimeter noise level,  number and length of 1/f baselines and total number of samples.  It can be obtained by using function build_noise_properties.\nrank: the rank of the current MPI process\ncomm: the MPI communicator object.\n\nThe following arguments are optional:\n\nthreshold is used by the conjugated-gradient algorithm. When the  residual error of the iteration goes below this value, the iteration  stops. The smaller the value, the more accurate the solution.\nmax_iter is the maximum number of iterations to be executed in the  conjugated-gradient algorithm. If the algorithm does not converge  within this number of iterations, the process will quit without  having reached the convergence threshold (see the threshold  keyword above).\nIf save_baseline_history is true, the return value will contain the sequence of baselines tested by the CG algorithm. Each MPI process will hold its own baselines.\n\nThe function returns a DestripingResults object containings the destriped map, the sequence of baselines, and other information describing the convergence of the CG algorithm.\n\nRemarks\n\nThe length of the arrays pix_idx and tod must be the same;\nIf you do not specify comm, no MPI will be used\n\n\n\n\n\n"
},

{
    "location": "mapmaking/#Map-making-functions-1",
    "page": "Map making",
    "title": "Map-making functions",
    "category": "section",
    "text": "Stripeline implements a few functions to produce maps out of time-ordered data.tod2map_mpi\nbaseline2map_mpi\ndestripe"
},

{
    "location": "simulation_tutorial/#",
    "page": "Simulation tutorial",
    "title": "Simulation tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "simulation_tutorial/#Simulation-tutorial-1",
    "page": "Simulation tutorial",
    "title": "Simulation tutorial",
    "category": "section",
    "text": "This aim of this tutorial is to describe how you can use the functions of Stripeline repository to perform a complete simulation of the LSPE/STRIP experiment.So far, the simulation includes:simulation of the telescope scanning the sky\nsimulation of noise (white and 1/f)\nproduction of a tod \nproduction of a map by using the destriping techniqueTwo examples will be presented: a simple case, in which we produce and analyze a very small TOD (1 day observation).  You should refer to this case if you want to use a Jupyter Notebook to perform preliminary studies (since you cannot use MPI in a notebook).\na more general and realistic case, suitable also for the production of large TODs.  In this case it will be necessary to use MPI functions."
},

{
    "location": "simulation_tutorial/#.-Simple-Case:-small-TOD-1",
    "page": "Simulation tutorial",
    "title": "1. Simple Case: small TOD",
    "category": "section",
    "text": "First of all, you should import the following packages:import Healpix\nimport Random\nimport CorrNoise\nimport Stripeline\nconst Sl = Stripeline\nusing FITSIOIn this simple case we can avoid using MPI, but we need the following line:comm = missing\nThis is because some of the functions we will use (e.g., the destripe function) require the MPI communicator as input parameter.Then, you can start by setting the simulation parameters:the number of days of observation\nthe number of polarimeters to simulate\nthe sampling frequency (Hz)\nthe length of 1/f baselines (s) as input for the destriper\nthe NSIDE you prefer for your output map\nthe temperature of the sky signal (K)\nthe temperature of the atmosphere (K)\nthe temperature of the telescope (K)\nthe noise temperature of the polarimeters (K)\nthe knee frequency of the polarimeters (Hz)\nthe bandwidth of the polarimeters (Hz)N.B. In the current version of the simulation we consider all polarimeters identical in properties.#Simulation parameters\n\ndays_of_observation = 1\nnum_of_polarimeters = 1\nbaseline_length_s = 10  \nfsamp_hz = 40\nNSIDE = 256\n\ntcmb_k = 3\ntatm_k = 15\nttel_k = 3\ntnoise_k = 35\nfknee_hz = 0.01\nβ_hz = 7e9\nAt this point, you can compute some parameters we will need later: the total observation time (in s), the integration time, the total system temperature and the receiver sensitivity.total_time = days_of_observation * 24 * 3600\ntsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k\nτ_s = 1 / fsamp_hz\nσ_k = (tsys_k / sqrt(β_hz * τ_s))Open now the input map, that is to say the sky you want to scan with your instrument.  You can find two example input maps (with two different resolutions) here.Those maps have been produced with PySM and are specific for the STRIP case: they are 43 GHz maps of polarized emission only (cmb, synchrotron and dust).  If you want to produce your own input map, you can use this python script.inputmap = Healpix.readMapFromFITS(raw\"PySM_inputmap_nside256.fits\", 2 , Float64)\ninputmap_resol = inputmap.resolution\n\nresol = Healpix.Resolution(NSIDE) #desired resolution for output map\nnum_of_pixels = resol.numOfPixelsYou can now scan the input map according to your scanning strategy, thus producing the pure signal \"sky TOD\".  In this case, we use the nominal scanning strategy for STRIP (Tenerife latitude, 20° of elevation angle, 1 RPM)#Generate sky tod\n\ntimes = 0:τ_s: (total_time-τ_s)\n(dirs, ψ) = Sl.genpointings([0., 0., 1.], times; latitude_deg=28.29) do time_s\n    return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1.))\nend\n\npix_idx_inputmap = Healpix.ang2pixRing.(Ref(inputmap_resol), dirs[:, 1], dirs[:, 2])\npix_idx = Healpix.ang2pixRing.(Ref(resol), dirs[:, 1], dirs[:, 2])\n\nsky_tod = inputmap.pixels[pix_idx_inputmap]Now you should add noise to your sky tod. We simulate both white noise and 1/f noise, in accordance with the noise properties of the polarimeters specified at the beginning of the script.  To do that, we use the functions of the module CorrNoise.jl.#Generate noise\n\nseed = rand(1:1000)\nrng = CorrNoise.OofRNG(CorrNoise.GaussRNG(Random.MersenneTwister(seed)), -1, 1.15e-5, fknee_hz, fsamp_hz);\nnoise_tod = [CorrNoise.randoof(rng) * σ_k for i in 1:(fsamp_hz * total_time)]\nFinally, you can get the final, realistic TOD just by doing:tod = sky_tod + noise_todOnce you have completed the simulation, you can do data analysis! For instance, you can call the destriper and clean the map from 1/f noise.(N.B. the destriper needs in input an array containing the lengths of all 1/f baselines.  For the sake of semplicity, we consider baselines of equal length).#Run the destriper\n\nbaseline_len = repeat([baseline_length_s*fsamp_hz],Int64(total_time/baseline_length_s))\n\n(destr_map, a) = Sl.destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)The output of the destriper are:destr_map : the destriped map, cleaned from 1/f noise.\na : the baselines array.If you wish, you can finally convert the destriped map (Array{Float64,1}) into a HEALPix map and save it in a FITS file:#save file \n\nmapfile = Healpix.Map{Float64,Healpix.RingOrder}(NSIDE)\nmapfile.pixels = destr_map\n\nHealpix.saveToFITS(mapfile, \"destriped_map.fits\", typechar = \"D\")To run this script you can use the following command:julia simplecase.jl"
},

{
    "location": "simulation_tutorial/#.-General-Case-1",
    "page": "Simulation tutorial",
    "title": "2. General Case",
    "category": "section",
    "text": "In realistic cases, TODs are really huge (millions or billions of samples!).  This means a lot of memory allocation. In the case of Strip, 49 polarimeters observing the sky for 2 years with a sampling frequency of 100 Hz produce 300 billions of samples, which means about 2 Terabytes of memory allocation!It is much more than a single computer can support. It is thus compulsory to split the TOD simulation and analysis between different computing units, by using MPI.Let\'s go through all the points of the simulation and see how things change:You have to add the MPI package to your dependencies. import MPI\n\nimport Healpix\nimport Random\nimport CorrNoise\nimport Stripeline\nconst Sl = Stripeline\nusing FITSIO\nNothing changes with respect to the case above.#Simulation parameters\n\ndays_of_observation = 1\nnum_of_polarimeters = 1\nbaseline_length_s = 10  \nfsamp_hz = 40\nNSIDE = 256\n\ntcmb_k = 3\ntatm_k = 15\nttel_k = 3\ntnoise_k = 35\nfknee_hz = 0.01\nβ_hz = 7e9\nYou should add the computation of the total number of samples per polarimeter and total number of baselines for polarimeter, which we will need later.total_time = days_of_observation * 24 * 3600\ntsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k\nτ_s = 1 / fsamp_hz\nσ_k = (tsys_k / sqrt(β_hz * τ_s))\n\nsamples_per_pol = total_time*fsamp_hz \nbaselines_per_pol = Int64(total_time/baseline_length_s)Nothing changes.inputmap = Healpix.readMapFromFITS(raw\"PySM_inputmap_nside256.fits\", 2 , Float64)\ninputmap_resol = inputmap.resolution\n\nresol = Healpix.Resolution(NSIDE) #desired resolution for output map\nnum_of_pixels = resol.numOfPixelsBefore scanning the input map we have to conveniently split the TOD production among the available computing units.  First of all, we need to initialize MPI:MPI.Init()\n\ncomm = MPI.COMM_WORLD\nrank = MPI.Comm_rank(comm)\ncommsize = MPI.Comm_size(comm) \nThen, we can proceed with the TOD splitting.The splitting is done in a way that each rank gets a whole number of 1/f baselines to simulate and that the TOD chunks are of as similar length as possible.We also need to tell to each unit which detector to simulate and from which to which time.  By using the get_chunk_properties function, we can obtain these information for the current rank.#Split tod production \n\nbaselines_per_process = Sl.split_into_n(num_of_polarimeters*baselines_per_pol, commsize)\nchunks = Sl.split_tod_mpi(total_time, baseline_length_s, baselines_per_process, commsize)\nthis_rank_chunk = chunks[rank+1]\n\n(detector_number, first_time, last_time, num_of_baselines, num_of_samples) = Sl.get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)\nConcerning the production of the sky TOD, of course each computing unit will produce its own partial TOD using the information obtained before.  A loop on detectors has been added, since in the general case the simulation involves more than one detector, and moreover, each rank may have to simulate partial tods for different detectors.#Generate sky tod\n\npix_idx = Int64[]\npix_idx_inputmap = Int64[]\nsky_tod = Float64[]\n\n\nfor i in 1:length(this_rank_chunk)   #loop on detectors\n\n    #generate pointings according to STRIP scanning strategy\n\n    times = first_time[i]:τ_s:last_time[i]\n\n \n    (dirs, ψ) = Sl.genpointings([0., 0., 1.], times; latitude_deg=28.29) do time_s\n        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1.))\n    end\n\n    partial_pix_idx_inputmap = Healpix.ang2pixRing.(Ref(inputmap_resol), dirs[:, 1], dirs[:, 2])\n    partial_pix_idx = Healpix.ang2pixRing.(Ref(resol), dirs[:, 1], dirs[:, 2])\n\n    #build the sky tod\n    partial_sky_tod = inputmap.pixels[partial_pix_idx_inputmap]\n    global sky_tod = append!(sky_tod, partial_sky_tod)\n    global pix_idx = append!(pix_idx, partial_pix_idx)\nend\nTo generate noise, you e can use the generate_noise_mpi function, which directly returns the partial noise TOD for the current rank.#Generate noise\n\nnoise_tod = Sl.generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, fsamp_hz, σ_k, fknee_hz, rank, comm)nothing changes.tod = sky_tod + noise_todnothing changes apart from baseline_len definition (since now different units can have a different amount of baselines to compute).The destripe function already takes into account the presence of multiple computing units: each partial TOD is loaded separately and partial maps are build, but then MPI functions are called to make different ranks \"talk together\" in order to obtain in output a single global destriped map.#Run the destriper\n\nbaseline_len = repeat([baseline_length_s*fsamp_hz], baselines_per_process[rank+1])\n(destr_map, a) = Sl.destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)If you want to save the destriped map in a .fits file, you should make just one rank do that.#save file \n\nif(rank==0)\n    mapfile = Healpix.Map{Float64,Healpix.RingOrder}(NSIDE)\n    mapfile.pixels = destr_map\n    Healpix.saveToFITS(mapfile, \"destriped_map.fits\", typechar = \"D\")\nendFinally, end your script by terminating the calling to MPI environment.MPI.Finalize()To run this script, you can use the following code (using 3 MPI processes):mpirun -n 3 julia generalcase.jl"
},

]}
