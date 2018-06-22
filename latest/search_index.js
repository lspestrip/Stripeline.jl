var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "index.html#Stripeline-User\'s-Manual-1",
    "page": "Introduction",
    "title": "Stripeline User\'s Manual",
    "category": "section",
    "text": "An implementation of a simulation/data analysis pipeline for the LSPE/STRIP instrument.To install Stripeline, start Julia and type the following command:Pkg.clone(\"https://github.com/lspestrip/Stripeline\")To run the test suite, type the following command:Pkg.test(\"Stripeline\")In this manual, we will often assume that Stripeline has been imported using the following commands:import Stripeline\nconst Sl = StripelineIn this way, we can call functions like genpointings using the syntax Sl.genpointings, instead of the longer Stripeline.genpointings."
},

{
    "location": "index.html#Documentation-1",
    "page": "Introduction",
    "title": "Documentation",
    "category": "section",
    "text": "The documentation was built using Documenter.jl.println(\"Documentation built $(now()) with Julia $(VERSION).\") # hide"
},

{
    "location": "index.html#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "instrumentdb.html#",
    "page": "Instrument database",
    "title": "Instrument database",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "instrumentdb.html#Instrument-database-1",
    "page": "Instrument database",
    "title": "Instrument database",
    "category": "section",
    "text": "InstrumentDB takes advantage of the structures Detector and Horn to retrieve information about feed horns and detectors from a YAML file. There are a set of YAML files containing the default configuration for the STRIP instrument in the repository."
},

{
    "location": "instrumentdb.html#Quick-introduction-1",
    "page": "Instrument database",
    "title": "Quick introduction",
    "category": "section",
    "text": "The following example initializes an object of type InstrumentDB with the values referred to the standard STRIP instrument:using Stripeline; # hide\ndb = InstrumentDB()As db is a struct, its field can be accessed with the usual dot notation. The two fields in db are focalplane and detectors. They are both dictionaries, associating horn names to Horn objects and detectors IDs to Detector objects, respectively:db.focalplane[\"I0\"]\ndb.detectors[2]The structure Detector is complex, as it is built over three other structures:BandshapeInfo\nSpectrumInfo\nNoiseTemperatureInfoAll these structures know how to show themselves on the REPL:db.detectors[2].bandshape\ndb.detectors[2].spectrum\ndb.detectors[2].tnoiseFor more information about the fields in the structures listed above, as well as their meaning, keep reading."
},

{
    "location": "instrumentdb.html#Stripeline.InstrumentDB",
    "page": "Instrument database",
    "title": "Stripeline.InstrumentDB",
    "category": "type",
    "text": "STRIP instrument database\n\nThe \"database\" contains information about feed horns and polarimeters:\n\nThe field focalplane is a dictionary (mapping) associating the string identifying a horn (e.g., I0) with a Horn structure;\nThe field detectors is a dictionary associating the ID of the polarimeter (e.g., 2 stands for STRIP02) with a Detector structure.\n\nYou should usually create an object of this kind using the default constructor, which parses a set of YAML files containing the real parameters of the instrument.\n\nExamples\n\njulia> db = InstrumentDB();\n\njulia> print(\"Number of horns in the database: $(length(keys(db.focalplane)))\")\nNumber of horns in the database: 49\n\njulia> print(\"Number of polarimeters in the database: $(length(keys(db.detectors)))\")\nNumber of polarimeters in the database: 66\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.Horn",
    "page": "Instrument database",
    "title": "Stripeline.Horn",
    "category": "type",
    "text": "Information about a STRIP horn\n\nThis structure holds a number of parameters relative to each feed horn in the STRIP focal plane.\n\nYou should initialize Horn objects via the InstrumentDB constructor, which loads their definition from a STRIP instrument database in YAML format.\n\nField Type Meaning\nname String Name of the horn, e.g., I0\nid Int Unique number of the horn, starting from 1\npolid Int Unique ID of the polarimeter associated with the horn\nmoduleid Int Number of the horn within the module, from 0 to 6\ncolor String Name of the color associated with the module\norientation Array{Float64} 3D vector containing the orientation of the horn in the sky\nfwhm_x_deg Float64 FWHM of the beam along the X axis, in degrees\nfwhm_y_deg Float64 FWHM of the beam along the Y axis, in degrees\nmain_spillover Float64 Main reflector spillover\nsub_spillover Float64 Sub-reflector spillover\nxpd_db Float64 Cross-polarization, in dB\ndirectivity_dbi Float64 Directivity, in dBi\nellipticity Float64 Ellipticity\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.Detector",
    "page": "Instrument database",
    "title": "Stripeline.Detector",
    "category": "type",
    "text": "Information about a STRIP detector\n\nThis structure holds information about a STRIP polarimeter.\n\nYou should initialize Detector objects via the InstrumentDB constructor, which loads their definition from a local STRIP instrument database.\n\nField Type Meaning\nid Int Integer ID of the polarimeter, e.g., 2 for STRIP02\nname String Full name of the polarimeter, e.g., STRIP02\nband String Band: it can either be Q or W\nbandshape BandshapeInfo Information about the bandpass response\nspectrum SpectrumInfo Information about the noise spectrum (white noise and 1/f noise)\ntnoise NoiseTemperatureInfo Information about the noise temperature\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.BandshapeInfo",
    "page": "Instrument database",
    "title": "Stripeline.BandshapeInfo",
    "category": "type",
    "text": "BandshapeInfo\n\nInformation about the spectral band response of a polarimeter.\n\nField Type Meaning\ncenter_frequency_hz Float64 Estimate for the center frequency, in Hz\ncenter_frequency_err_hz Float64 Estimated error on the center frequency, in Hz\nbandwidth_hz Float64 Estimated bandwidth, in Hz\nbandwidth_err_hz Float64 Estimated error on the bandwidth, in Hz\nlowest_frequency_hz Float64 Lowest frequency of the bandshape in response, in Hz\nhighest_frequency_hz Float64 Highest frequency of the bandshape in response, in Hz\nnum_of_frequencies Int Number of samples in response\nresponse Array{1,Float64} Profile of the bandshape (pure numbers)\ntest_id Int ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.SpectrumInfo",
    "page": "Instrument database",
    "title": "Stripeline.SpectrumInfo",
    "category": "type",
    "text": "SpectrumInfo\n\nInformation about the noise spectrum of the output of a polarimeter.\n\nField Type Meaning\nslope_i Float64 The slope (alpha) of the 1/f component of the noise in the I signal\nslope_i_err Float64 Error associated with the value of slope_i\nslope_q Float64 Same as slope_i, but for the Q signal\nslope_q_err Float64 Error associated with the value of slope_q\nslope_u Float64 Same as slope_i, but for the U signal\nslope_u_err Float64 Error associated with the value of slope_u\nfknee_i_hz Float64 Knee frequency of the I signal, in Hz\nfknee_i_err_hz Float64 Error associated with the value of fknee_i_hz\nfknee_q_hz Float64 Knee frequency of the Q signal, in Hz\nfknee_q_err_hz Float64 Error associated with the value of fknee_q_hz\nfknee_u_hz Float64 Knee frequency of the U signal, in Hz\nfknee_u_err_hz Float64 Error associated with the value of fknee_u_hz\nwn_i_k2_hz Float64 White noise level for the I signal, in K^2 Hz\nwn_i_err_k2_hz Float64 Error associated with the value of wn_i_k2_hz\nwn_q_k2_hz Float64 White noise level for the Q signal, in K^2 Hz\nwn_q_err_k2_hz Float64 Error associated with the value of wn_q_k2_hz\nwn_u_k2_hz Float64 White noise level for the U signal, in K^2 Hz\nwn_u_err_k2_hz Float64 Error associated with the value of wn_u_k2_hz\ntest_id Int ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.NoiseTemperatureInfo",
    "page": "Instrument database",
    "title": "Stripeline.NoiseTemperatureInfo",
    "category": "type",
    "text": "NoiseTemperatureInfo\n\nInformation about the noise temperature of a polarimeter. This structure is used for the field tnoise of the Detector struct.\n\nField Type Meaning\ntnoise_k Float64 Noise temperature computed from tnoise_values_k, in K\ntnoise_err_k Float64 Error associated with tnoise_k, computed from tnoise_values_k\ntest_ids Array{Int,1} List of unit-level test IDs used to estimate the noise temperature\nanalysis_ids Array{Int,1} List of unit-level analysis report IDs used to estimate the noise temperature\nvalues_k Array{Float64,1} List of noise temperatures estimated from the tests\n\n\n\n"
},

{
    "location": "instrumentdb.html#Structures-1",
    "page": "Instrument database",
    "title": "Structures",
    "category": "section",
    "text": "InstrumentDB\nHorn\nDetector\nBandshapeInfo\nSpectrumInfo\nNoiseTemperatureInfo"
},

{
    "location": "instrumentdb.html#Stripeline.defaultdbfolder",
    "page": "Instrument database",
    "title": "Stripeline.defaultdbfolder",
    "category": "function",
    "text": "defaultdbfolder()\n\nReturn a string containing the (local) full path to the YAML files containing the reference instrument DB.\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.parsefpdict",
    "page": "Instrument database",
    "title": "Stripeline.parsefpdict",
    "category": "function",
    "text": "parsefpdict(fpdict)\n\nReturn a dictionary associating an horn name (e.g., I0) to a Horn object containing information about some horn in the STRIP focal plane. The information are parsed from fpdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_focal_plane.yaml\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.parsedetdict",
    "page": "Instrument database",
    "title": "Stripeline.parsedetdict",
    "category": "function",
    "text": "parsedetdict(detdict)\n\nReturn a dictionary associating an integer number to a Detector object containing information about the STRIP detector with the corresponding number. The information are parsed from detdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_detectors.yaml\n\n\n\n"
},

{
    "location": "instrumentdb.html#Loading-custom-databases-1",
    "page": "Instrument database",
    "title": "Loading custom databases",
    "category": "section",
    "text": "It is not needed to load the default instrument database, as Stripeline provides a number of additional functions to build mock databases from dictionaries.defaultdbfolder\nparsefpdict\nparsedetdict"
},

{
    "location": "rng.html#",
    "page": "Random numbers",
    "title": "Random numbers",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "rng.html#Random-number-generators-1",
    "page": "Random numbers",
    "title": "Random number generators",
    "category": "section",
    "text": "Although Julia already implements a number of pseudo-random number generators, Stripeline implements its own generators. Strip re-implements the same generators used in the pipeline of the Planck/LFI instrument, which provided several types of distributions:Uniform distribution (Flat128RNG), with period 2^128;\nGaussian distribution (GaussRNG);\n1f^2\ndistribution (Oof2RNG);\n1f^\ndistribution, with   2  (OofRNG).Each generator but Flat128RNG uses a simpler generator internally. This generator must sample from a given distribution, but it does not need to be a generator provided by Stripeline. For instance, the Gaussian generator GaussRNG employs an uniform generator, which can either be Flat128RNG or one of the generators provided by Julia like MersenneTwister. For instance, here is an example which shows how to use Flat128RNG:using Stripeline # hide\ngauss1 = GaussRNG(initflatrng128(1234))\nprint([randn(gauss1) for i in 1:4])We use initflatrng128, as it creates a Flat128RNG object with some sensible defaults (specifically, it is configured to produce the same sequence of random numbers as the ones produced by the Planck/HFI pipeline, if the seeds are the same). And here is the same example, using a MersenneTwister generator:gauss2 = GaussRNG(MersenneTwister(1234))\nprint([randn(gauss2) for i in 1:4])Of course, the numbers are different. They are however drawn from the same distribution (Gaussian curve with mean 0 and σ=1)."
},

{
    "location": "rng.html#Stripeline.Flat128RNG",
    "page": "Random numbers",
    "title": "Stripeline.Flat128RNG",
    "category": "type",
    "text": "Flat128RNG\n\nState of the base-128 uniform random generator. Initialize this using the function initflatrng128.\n\n\n\n"
},

{
    "location": "rng.html#Stripeline.initflatrng128",
    "page": "Random numbers",
    "title": "Stripeline.initflatrng128",
    "category": "function",
    "text": "initflatrng128(xstart = 123456789, ystart = 362436069, zstart = 521288629, wstart = 88675123)\n\nInitialize a flat random number generator with period 2^128. To draw random numbers, use the Base.rand function as usual.\n\nExample:\n\nrng = initflatrng128()\nprint([rand(rng) for i in 1:4])\n\n\n\n"
},

{
    "location": "rng.html#Uniform-generator-1",
    "page": "Random numbers",
    "title": "Uniform generator",
    "category": "section",
    "text": "Flat128RNG\ninitflatrng128"
},

{
    "location": "rng.html#Stripeline.GaussRNG",
    "page": "Random numbers",
    "title": "Stripeline.GaussRNG",
    "category": "type",
    "text": "GaussRNG(flatrng::AbstractRNG)\nGaussRNG(seed=0)\n\nInitialize a Gaussian RNG. The parameter flatrng must be a uniform RNG. If a seed is used, then a MersenneTwister RNG is used.\n\n\n\n"
},

{
    "location": "rng.html#Gaussian-generator-1",
    "page": "Random numbers",
    "title": "Gaussian generator",
    "category": "section",
    "text": "GaussRNG"
},

{
    "location": "rng.html#Stripeline.Oof2RNG",
    "page": "Random numbers",
    "title": "Stripeline.Oof2RNG",
    "category": "type",
    "text": "Oof2RNG(normrng, fmin, fknee, fsample)\n\nCreate a Oof2RNG RNG object. It requires a gaussian RNG generator in normrng (use GaussRNG), the minimum frequency (longest period) in fmin, the knee frequency and the sampling frequency. The measure unit of the three frequencies must be the same (e.g., Hz).\n\nUse randoof2 to draw samples from a Oof2RNG object, as in the following example:\n\nrng = Oof2RNG(GaussRNG(), 1e-3, 1.0, 1e2)\nprint([randoof2(rng) for i in 1:4])\n\n\n\n"
},

{
    "location": "rng.html#1/f2-generator-1",
    "page": "Random numbers",
    "title": "1f^2 generator",
    "category": "section",
    "text": "Oof2RNG"
},

{
    "location": "rng.html#Stripeline.OofRNG",
    "page": "Random numbers",
    "title": "Stripeline.OofRNG",
    "category": "type",
    "text": "OofRNG(normrng, fmin, fknee, fsample)\n\nCreate a OofRNG RNG object. It requires a gaussian RNG generator in normrng (use GaussRNG), the slope α of the noise in slope, the minimum frequency (longest period) in fmin, the knee frequency and the sampling frequency. The measure unit of the three frequencies must be the same (e.g., Hz).\n\nUse randoof to draw samples from a OofRNG object, as in the following example:\n\nrng = OofRNG(GaussRNG(), 1.5, 1e-3, 1.0, 1e2)\nprint([randoof(rng) for i in 1:4])\n\n\n\n"
},

{
    "location": "rng.html#1/fα-generator-(with-α-2)-1",
    "page": "Random numbers",
    "title": "1f^ generator (with   2)",
    "category": "section",
    "text": "OofRNG"
},

{
    "location": "scanning.html#",
    "page": "Scanning strategy",
    "title": "Scanning strategy",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "scanning.html#Stripeline.genpointings",
    "page": "Scanning strategy",
    "title": "Stripeline.genpointings",
    "category": "function",
    "text": "genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)\n\nGenerate a set of pointings for some STRIP detector. The parameter wheelanglesfn must be a function which takes as input a time in seconds and returns a 3-tuple containing the angles (in radians) of the three motors:\n\nThe boresight motor\nThe altitude motor\nThe ground motor\n\nThe parameter dir must be a normalized vector which tells the pointing direction of the beam (boresight is [0, 0, 1]). The parameter timerange_s is either a range or a vector which specifies at which times (in second) the pointings should be computed. The keyword latitude_deg should contain the latitude (in degrees, N is positive) of the location where the observation is made.\n\nReturn a 2-tuple containing the directions (a N×2 array containing the colatitude and the longitude) and the polarization angles at each time step. Directions are expressed in equatorial coordinates.\n\nExample:\n\ngenpointings([0, 0, 1], 0:0.1:1) do time_s\n    # Boresight motor keeps a constant angle equal to 0°\n    # Altitude motor remains at 20° from the Zenith\n    # Ground motor spins at 1 RPM\n    return (0.0, deg2rad(20.0), timetorotang(time_s, 1))\nend\n\n\n\n"
},

{
    "location": "scanning.html#Stripeline.timetorotang",
    "page": "Scanning strategy",
    "title": "Stripeline.timetorotang",
    "category": "function",
    "text": "timetorotang(time, rpm)\n\nConvert a time into a rotation angle, given the number of rotations per minute. The time should be expressed in seconds. The return value is in radians. time can either be a scalar or a vector.\n\n\n\n"
},

{
    "location": "scanning.html#Simulating-the-scanning-strategy-1",
    "page": "Scanning strategy",
    "title": "Simulating the scanning strategy",
    "category": "section",
    "text": "genpointings\ntimetorotang"
},

]}
