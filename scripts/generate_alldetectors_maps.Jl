using Stripeline, Healpix, ColorSchemes, Plots, FITSIO
const Sl = Stripeline
const Hp = Healpix

pars_arg = ARGS

elevation_angle_deg = parse(Float64, pars_arg[1])

g = FITS(join([pars_arg[2], "I0_", Int(elevation_angle_deg), "°.fits"]), "r")
h = read_header(g[2])

hits_map = Healpix.Map{Float64, Healpix.RingOrder}(h["NSIDE"])
hits_map.pixels = zeros(length(hits_map.pixels))

db = Sl.InstrumentDB()
for fp in db.focalplane

    load_path = join([pars_arg[2],
                      fp[2].name,
                      "_",
                      Int(elevation_angle_deg),
                      "°.fits"])

    single_detector_map = Healpix.readMapFromFITS(load_path, 1, Float64)
    
    hits_map.pixels += single_detector_map.pixels

end

save_to_path = join([pars_arg[2],
                     "alldetectors_",
                     Int(elevation_angle_deg),
                      "°.fits"])

Healpix.saveToFITS(hits_map, save_to_path, typechar = "D")

f = FITS(save_to_path, "r+")

write_key(f[2], "LATITUDE", h["LATITUDE"], "Tenerife Latitude [Deg]")
write_key(f[2], "LONGITUD", h["LONGITUD"], "Tenerife Longitude [Deg]")
write_key(f[2], "HEIGHT", h["HEIGHT"], "Tenerife Height [m]")

write_key(f[2], "DURATION", h["DURATION"], "[sec]")
write_key(f[2], "TSTART", h["TSTART"],
          "Starting Date and time of the observations")    
write_key(f[2], "TSTOP", h["TSTOP"], "Ending Date and time of the observations")

write_key(f[2], "ANG_EL", h["ANG_EL"], "Elevation angle [Deg]")
write_key(f[2], "SPIN_VEL", h["SPIN_VEL"], "Spin velocity [rpm]")
write_key(f[2], "S_RATE", h["S_RATE"], "Sampling rate [Hz]")

write_key(f[2], "DETECTOR", "all detectors", "Sum of all detectors")

close(f)

