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
    "text": "STRIP instrument database\n\nThe \"database\" contains information about feed horns and polarimeters:\n\nThe field focalplane is a dictionary (mapping) associating the string identifying a horn (e.g., I0) with a Horn structure;\nThe field detectors is a dictionary associating the ID of the polarimeter (e.g., 2 stands for STRIP02) with a Detector structure.\n\nYou should usually create an object of this kind using the default constructor, which parses a set of YAML files containing the real parameters of the instrument.\n\nExamples\n\njulia> db = InstrumentDB();\n\njulia> print(\"Number of horns in the database: $(length(keys(db.focalplane)))\")\nNumber of horns in the database: 49\n\njulia> print(\"Number of polarimeters in the database: $(length(keys(db.detectors)))\")\nNumber of polarimeters in the database: 66\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.Horn",
    "page": "Instrument database",
    "title": "Stripeline.Horn",
    "category": "type",
    "text": "Information about a STRIP horn\n\nThis structure holds a number of parameters relative to each feed horn in the STRIP focal plane.\n\nYou should initialize Horn objects via the InstrumentDB constructor, which loads their definition from a STRIP instrument database in YAML format.\n\nField Type Meaning\nname String Name of the horn, e.g., I0\nid Int Unique number of the horn, starting from 1\npolid Int Unique ID of the polarimeter associated with the horn\nmoduleid Int Number of the horn within the module, from 0 to 6\ncolor String Name of the color associated with the module\norientation Array{Float64} 3D vector containing the orientation of the horn in the sky\nfwhm_x_deg Float64 FWHM of the beam along the X axis, in degrees\nfwhm_y_deg Float64 FWHM of the beam along the Y axis, in degrees\nmain_spillover Float64 Main reflector spillover\nsub_spillover Float64 Sub-reflector spillover\nxpd_db Float64 Cross-polarization, in dB\ndirectivity_dbi Float64 Directivity, in dBi\nellipticity Float64 Ellipticity\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.Detector",
    "page": "Instrument database",
    "title": "Stripeline.Detector",
    "category": "type",
    "text": "Information about a STRIP detector\n\nThis structure holds information about a STRIP polarimeter.\n\nYou should initialize Detector objects via the InstrumentDB constructor, which loads their definition from a local STRIP instrument database.\n\nField Type Meaning\nid Int Integer ID of the polarimeter, e.g., 2 for STRIP02\nname String Full name of the polarimeter, e.g., STRIP02\nband String Band: it can either be Q or W\nbandshape BandshapeInfo Information about the bandpass response\nspectrum SpectrumInfo Information about the noise spectrum (white noise and 1/f noise)\ntnoise NoiseTemperatureInfo Information about the noise temperature\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.BandshapeInfo",
    "page": "Instrument database",
    "title": "Stripeline.BandshapeInfo",
    "category": "type",
    "text": "BandshapeInfo\n\nInformation about the spectral band response of a polarimeter.\n\nField Type Meaning\ncenter_frequency_hz Float64 Estimate for the center frequency, in Hz\ncenter_frequency_err_hz Float64 Estimated error on the center frequency, in Hz\nbandwidth_hz Float64 Estimated bandwidth, in Hz\nbandwidth_err_hz Float64 Estimated error on the bandwidth, in Hz\nlowest_frequency_hz Float64 Lowest frequency of the bandshape in response, in Hz\nhighest_frequency_hz Float64 Highest frequency of the bandshape in response, in Hz\nnum_of_frequencies Int Number of samples in response\nresponse Array{1,Float64} Profile of the bandshape (pure numbers)\ntest_id Int ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.SpectrumInfo",
    "page": "Instrument database",
    "title": "Stripeline.SpectrumInfo",
    "category": "type",
    "text": "SpectrumInfo\n\nInformation about the noise spectrum of the output of a polarimeter.\n\nField Type Meaning\nslope_i Float64 The slope (alpha) of the 1/f component of the noise in the I signal\nslope_i_err Float64 Error associated with the value of slope_i\nslope_q Float64 Same as slope_i, but for the Q signal\nslope_q_err Float64 Error associated with the value of slope_q\nslope_u Float64 Same as slope_i, but for the U signal\nslope_u_err Float64 Error associated with the value of slope_u\nfknee_i_hz Float64 Knee frequency of the I signal, in Hz\nfknee_i_err_hz Float64 Error associated with the value of fknee_i_hz\nfknee_q_hz Float64 Knee frequency of the Q signal, in Hz\nfknee_q_err_hz Float64 Error associated with the value of fknee_q_hz\nfknee_u_hz Float64 Knee frequency of the U signal, in Hz\nfknee_u_err_hz Float64 Error associated with the value of fknee_u_hz\nwn_i_k2_hz Float64 White noise level for the I signal, in K^2 Hz\nwn_i_err_k2_hz Float64 Error associated with the value of wn_i_k2_hz\nwn_q_k2_hz Float64 White noise level for the Q signal, in K^2 Hz\nwn_q_err_k2_hz Float64 Error associated with the value of wn_q_k2_hz\nwn_u_k2_hz Float64 White noise level for the U signal, in K^2 Hz\nwn_u_err_k2_hz Float64 Error associated with the value of wn_u_k2_hz\ntest_id Int ID of the unit-level test used to characterize the bandshape\nanalysis_id Int ID of the unit-level analysis used to characterize the bandshape\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.NoiseTemperatureInfo",
    "page": "Instrument database",
    "title": "Stripeline.NoiseTemperatureInfo",
    "category": "type",
    "text": "NoiseTemperatureInfo\n\nInformation about the noise temperature of a polarimeter. This structure is used for the field tnoise of the Detector struct.\n\nField Type Meaning\ntnoise_k Float64 Noise temperature computed from tnoise_values_k, in K\ntnoise_err_k Float64 Error associated with tnoise_k, computed from tnoise_values_k\ntest_ids Array{Int,1} List of unit-level test IDs used to estimate the noise temperature\nanalysis_ids Array{Int,1} List of unit-level analysis report IDs used to estimate the noise temperature\nvalues_k Array{Float64,1} List of noise temperatures estimated from the tests\n\n\n\n\n\n"
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
    "text": "defaultdbfolder()\n\nReturn a string containing the (local) full path to the YAML files containing the reference instrument DB.\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.parsefpdict",
    "page": "Instrument database",
    "title": "Stripeline.parsefpdict",
    "category": "function",
    "text": "parsefpdict(fpdict)\n\nReturn a dictionary associating an horn name (e.g., I0) to a Horn object containing information about some horn in the STRIP focal plane. The information are parsed from fpdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_focal_plane.yaml\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Stripeline.parsedetdict",
    "page": "Instrument database",
    "title": "Stripeline.parsedetdict",
    "category": "function",
    "text": "parsedetdict(detdict)\n\nReturn a dictionary associating an integer number to a Detector object containing information about the STRIP detector with the corresponding number. The information are parsed from detdict, which should be a dictionary loaded from a YAML file. The default YAML file to be used is located in the folder returned by defaultdbfolder and is usually named strip_detectors.yaml\n\n\n\n\n\n"
},

{
    "location": "instrumentdb.html#Loading-custom-databases-1",
    "page": "Instrument database",
    "title": "Loading custom databases",
    "category": "section",
    "text": "It is not needed to load the default instrument database, as Stripeline provides a number of additional functions to build mock databases from dictionaries.defaultdbfolder\nparsefpdict\nparsedetdict"
},

{
    "location": "scanning.html#",
    "page": "Scanning strategy",
    "title": "Scanning strategy",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Stripeline\nend"
},

{
    "location": "scanning.html#Stripeline.genskypointings",
    "page": "Scanning strategy",
    "title": "Stripeline.genskypointings",
    "category": "function",
    "text": "genskypointings(t_start, t_stop, dirs; latitude_deg=0.0, longitude_deg=0.0, height_m=0.0)\n\nProject a set of pointings in the sky.  \n\nThe parameter t_start and t_start must be two strings which tells the exact  UTC date and time of the observation in \"iso\" format. The parameter dirs must  be a N×2 array containing the observed directions expressed in colatitude and  the longitude. The keywords latitude_deg, longitude_deg and height_m should contain the latitude (in degrees, N is positive), the longitude (in degrees,  counterclockwise is positive) and the height (in meters)  of the location where  the observation is made.\n\nReturn a 2-tuple containing the observed directions projected in sky coordinates  (a N×2 array containing the Right ascension and the Declination, in radians) at  each time step.  Directions are expressed in ICRS coordinates.\n\nExample:\n\ngenskypointings(\"2019-01-01 00:00:00\", \"2020-01-25 14:52:10.05\", dirs)\n\n\n\n\n\n"
},

{
    "location": "scanning.html#Stripeline.genpointings",
    "page": "Scanning strategy",
    "title": "Stripeline.genpointings",
    "category": "function",
    "text": "genpointings(wheelanglesfn, dir, timerange_s; latitude_deg=0.0)\n\nGenerate a set of pointings for some STRIP detector. The parameter wheelanglesfn must be a function which takes as input a time in seconds and returns a 3-tuple containing the angles (in radians) of the three motors:\n\nThe boresight motor\nThe altitude motor\nThe ground motor\n\nThe parameter dir must be a normalized vector which tells the pointing direction of the beam (boresight is [0, 0, 1]). The parameter timerange_s is either a range or a vector which specifies at which times (in second) the pointings should be computed. The keyword latitude_deg should contain the latitude (in degrees, N is positive) of the location where the observation is made.\n\nReturn a 2-tuple containing the directions (a N×2 array containing the colatitude and the longitude) and the polarization angles at each time step. Directions are expressed in equatorial coordinates.\n\nExample:\n\ngenpointings([0, 0, 1], 0:0.1:1) do time_s\n    # Boresight motor keeps a constant angle equal to 0°\n    # Altitude motor remains at 20° from the Zenith\n    # Ground motor spins at 1 RPM\n    return (0.0, deg2rad(20.0), timetorotang(time_s, 1))\nend\n\n\n\n\n\n"
},

{
    "location": "scanning.html#Stripeline.timetorotang",
    "page": "Scanning strategy",
    "title": "Stripeline.timetorotang",
    "category": "function",
    "text": "timetorotang(time, rpm)\n\nConvert a time into a rotation angle, given the number of rotations per minute. The time should be expressed in seconds. The return value is in radians. time can either be a scalar or a vector.\n\n\n\n\n\n"
},

{
    "location": "scanning.html#Simulating-the-scanning-strategy-1",
    "page": "Scanning strategy",
    "title": "Simulating the scanning strategy",
    "category": "section",
    "text": "genskypointings\ngenpointings\ntimetorotang"
},

{
    "location": "simulation_tutorial.html#",
    "page": "Simulation tutorial",
    "title": "Simulation tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "simulation_tutorial.html#Simulation-tutorial-1",
    "page": "Simulation tutorial",
    "title": "Simulation tutorial",
    "category": "section",
    "text": "This aim of this tutorial is to describe how you can use the functions of Stripeline repository to perform a complete simulation of the LSPE/STRIP experiment. So far, the simulation includes:simulation of the telescope scanning the sky\nsimulation of noise (white and 1/f)\nproduction of a tod \nproduction of a map by using the destriping techniqueTwo examples will be presented: a simple case, in which we produce and analyze a very small TOD (1 day observation)\na more general and realistic case, suitable also for the production of large TODs.In this case it will be necessary to use MPI functions."
},

{
    "location": "simulation_tutorial.html#.-Simple-Case:-small-TOD-1",
    "page": "Simulation tutorial",
    "title": "1. Simple Case: small TOD",
    "category": "section",
    "text": "First of all, you should import the following packages:import Healpix\nimport Random\nimport CorrNoise\nimport Stripeline\nconst Sl = Stripeline\nusing FITSIO\nIn this simple case we can avoid using MPI, but we need the following line: comm = missing\nsince some of the functions we will use (e.g. the destripe function) require the MPI communicator as input parameter. Then, you can start by setting the simulation parameters:the number of days of observation\nthe number of polarimeters to simulate\nthe sampling frequency (Hz)\nthe length of 1/f baselines (s) as input for the destriper\nthe NSIDE you prefer for your output map\nthe temperature of the sky signal (K)\nthe temperature of the atmosphere (K)\nthe temperature of the telescope (K)\nthe noise temperature of the polarimeters (K)\nthe knee frequency of the polarimeters (Hz)\nthe bandwidth of the polarimeters (Hz)N.B. In the current version of the simulation we consider all polarimeters identical in properties.#Simulation parameters\n\ndays_of_observation = 1\nnum_of_polarimeters = 1\nbaseline_length_s = 10  \nfsamp_hz = 40\nNSIDE = 256\n\ntcmb_k = 3\ntatm_k = 15\nttel_k = 3\ntnoise_k = 35\nfknee_hz = 0.01\nβ_hz = 7e9\nAt this point, you can compute some parameters we will need later: the total observation time (in s), the integration time, the total system temperature and the receiver sensitivity. total_time = days_of_observation * 24 * 3600\ntsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k\nτ_s = 1 / fsamp_hz\nσ_k = (tsys_k / sqrt(β_hz * τ_s))Open now the input map, that is to say the sky you want to scan with your instrument. You can find two example input maps (with two different resolutions) here: https://github.com/lspestrip/Stripeline.jl/tree/master/test/testfiles.Those maps have been produced with PySM and are specific for the STRIP case: they are 43 GHz maps of polarized emission only (cmb, synchrotron and dust). If you want to produce your own input map, you can use this python script. inputmap = Healpix.readMapFromFITS(raw\"PySM_inputmap_nside256.fits\", 2 , Float64)\ninputmap_resol = inputmap.resolution\n\nresol = Healpix.Resolution(NSIDE) #desired resolution for output map\nnum_of_pixels = resol.numOfPixelsYou can now scan the input map according to your scanning strategy, thus producing the pure signal \"sky TOD\".In this case we use nominal scanning strategy for STRIP (Tenerife latitude, 20 degrees of elevation angle, 1 rpm)#Generate sky tod\n\ntimes = 0:τ_s: (total_time-τ_s)\n(dirs, ψ) = Sl.genpointings([0, 0, 1], times; latitude_deg=28.29) do time_s\n    return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1))\nend\n\npix_idx_inputmap = Healpix.ang2pixRing.(Ref(inputmap_resol), dirs[:, 1], dirs[:, 2])\npix_idx = Healpix.ang2pixRing.(Ref(resol), dirs[:, 1], dirs[:, 2])\n\nsky_tod = inputmap.pixels[pix_idx_inputmap]Now you should add noise to your sky tod. We simulate both white noise and 1/f noise, in accordance with the noise properties of the polarimeters specified at the beginning of the script.To do that, we use the functions of the module CorrNoise.jl.<<<<<<< HEAD\n#Generate noise\n=======\n#Generate Noise\n>>>>>>> bcb5ddeacb482f26f17d9f8b642719999d58f234\n\nseed = rand(1:1000)\nrng = CorrNoise.OofRNG(CorrNoise.GaussRNG(Random.MersenneTwister(seed)), -1, 1.15e-5, fknee_hz, fsamp_hz);\nnoise_tod = [CorrNoise.randoof(rng) * σ_k for i in 1:(fsamp_hz * total_time)]\nfinally, you can get the final, realistic TOD just by doing:tod = sky_tod + noise_todOnce simulated your data, you can now perform data analysis.You can call the destriper and clean the map from 1/f noise.(N.B. the destriper needs in input an array containing the lengths of all 1/f baselines. For the sake of semplicity, we consider baselines of equal length).#Run the destriper\n\nbaseline_len = repeat([baseline_length_s*fsamp_hz],Int64(total_time/baseline_length_s))\n\n(destr_map, a) = Sl.destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)The output of the destriper are:destr_map : the destriped map, cleaned from 1/f noise.\na : the baselines array.If you wish, you can finally save the destriped map in a .fits file: #save file \n\nmapfile = Healpix.Map{Float64,Healpix.RingOrder}(NSIDE)\nmapfile.pixels = destr_map\nHealpix.saveToFITS(mapfile, \"!results/destriped_map.fits\", typechar = \"D\")To run this script you can do:julia simplecase.jl"
},

{
    "location": "simulation_tutorial.html#.-General-Case-1",
    "page": "Simulation tutorial",
    "title": "2. General Case",
    "category": "section",
    "text": "In realistic cases, TODs are really huge (millions or billions of samples!). This means a lot of memory allocation.STRIP case:49 polarimeters * 2 years * 365 days * 86400 s * 100 Hz ≈ 300 billion sampleswhich means about 2 Terabytes of memory allocation!It is much more a single computer can support. It is thus compulsory to split the TOD simulation and analysis between different computing units, by using MPI.Let\'s go through all the points of the simulation and see how things change:You have to add the MPI package to your dependencies. import MPI\n\nimport Healpix\nimport Random\nimport CorrNoise\nimport Stripeline\nconst Sl = Stripeline\nusing FITSIO\nnothing changes.#Simulation parameters\n\ndays_of_observation = 1\nnum_of_polarimeters = 1\nbaseline_length_s = 10  \nfsamp_hz = 40\nNSIDE = 256\n\ntcmb_k = 3\ntatm_k = 15\nttel_k = 3\ntnoise_k = 35\nfknee_hz = 0.01\nβ_hz = 7e9\nyou should add the computation of the total number of samples per polarimeter and total number of baselines for polarimeter, which we will need later.samples_per_pol = total_time*fsamp_hz \nbaselines_per_pol = Int64(total_time/baseline_length_s)\n\ntotal_time = days_of_observation * 24 * 3600\ntsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k\nτ_s = 1 / fsamp_hz\nσ_k = (tsys_k / sqrt(β_hz * τ_s))nothing changes.inputmap = Healpix.readMapFromFITS(raw\"PySM_inputmap_nside256.fits\", 2 , Float64)\ninputmap_resol = inputmap.resolution\n\nresol = Healpix.Resolution(NSIDE) #desired resolution for output map\nnum_of_pixels = resol.numOfPixelsBefore scanning the input map we have to conveniently split the TOD production among the available computing units. First of all, we need to initialize MPI:MPI.Init()\n\ncomm = MPI.COMM_WORLD\nrank = MPI.Comm_rank(comm)\ncommsize = MPI.Comm_size(comm) \nThen, we can proceed with the TOD splitting.The splitting is done in a way that each rank gets a whole number of 1/f baselines to simulate and that the TOD chunks are of as similar length as possible.We also need to tell to each unit which detector to simulate and from which to which time. By using the get_chunk_properties function, we can obtain these information for the current rank. #Split tod production \n\nbaselines_per_process = Sl.split_into_n(num_of_polarimeters*baselines_per_pol, commsize)\nchunks = Sl.split_tod_mpi(total_time, baseline_length_s, baselines_per_process, commsize)\nthis_rank_chunk = chunks[rank+1]\n\n(detector_number, first_time, last_time, num_of_baselines, num_of_samples) = Sl.get_chunk_properties(chunks, baseline_length_s, fsamp_hz, rank)\nConcerning the production of the sky TOD, of course each computing unit will produce its own partial TOD using the information obtained before.A loop on detectors has been added, since in the general case the simulation involves more than one detector, and moreover, each rank may have to simulate partial tods for different detectors.#Generate sky tod\n\npix_idx = Int64[]\npix_idx_planck = Int64[]\nskytod_Q = Float64[]\n\n\nfor i in 1:length(this_rank_chunk)   #loop on detectors\n\n    #generate pointings according to STRIP scanning strategy\n\n    times = first_time[i]:τ_s:last_time[i]\n\n    ltimes= length(times)\n    (dirs, ψ) = Sl.genpointings([0, 0, 1], times; latitude_deg=28.29) do time_s\n        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1))\n    end\n\n    partial_pix_idx_CMB = Healpix.ang2pixRing.(Ref(CMB_map_resol), dirs[:, 1], dirs[:, 2])\n    partial_pix_idx = Healpix.ang2pixRing.(Ref(resol), dirs[:, 1], dirs[:, 2])\n\n    #build the sky tod\n    partial_skytod_Q = CMB_map.pixels[partial_pix_idx_CMB]\n    global skytod_Q = append!(skytod_Q, partial_skytod_Q)\n    global pix_idx = append!(pix_idx, partial_pix_idx)\nend\nTo generate noise, you e can use the generate_noise_mpi function, which directly returns the partial noise TOD for the current rank. #Generate noise\n\nnoise = Sl.generate_noise_mpi(chunks, baselines_per_process, baseline_length_s, fsamp_hz, σ_k, fknee_hz, rank, comm)nothing changes.tod = sky_tod + noise_todnothing changes apart from baseline_len definition (since now different units can have a different amount of baselines to compute).The destripe function already takes into account the presence of multiple computing units: each partial TOD is loaded separately and partial maps are build, but then MPI functions are called to make different ranks \"talk together\" in order to obtain in output a single global destriped map.#Run the destriper\n\nbaseline_len = repeat([baseline_length_s*fsamp_hz], baselines_per_process[rank+1])\n(destr_map, a) = Sl.destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)If you want to save the destriped map in a .fits file, you should make just one rank do that. #save file \n\nif(rank==0)\n    mapfile = Healpix.Map{Float64,Healpix.RingOrder}(NSIDE)\n    mapfile.pixels = destr_map\n    Healpix.saveToFITS(mapfile, \"!results/destriped_map.fits\", typechar = \"D\")\nendFinally, end your script by terminating the calling to MPI environment.MPI.Finalize()To run this script you can do (e.g. 3 computing units)mpirun -n 3 julia generalcase.jl"
},

]}
