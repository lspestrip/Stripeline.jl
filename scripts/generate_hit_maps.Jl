using Stripeline, Healpix, ColorSchemes, Plots, FITSIO, Dates, MPI
const Sl = Stripeline

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm) 

function get_hits_map(elevation_angle::Float64,
                      spin_velocity::Float64,
                      sampling_rate::Float64,
                      time_duration::Int,
                      t_start::DateTime,
                      t_stop::DateTime,
                      NSIDE::Int,
                      latitude::Float64,
                      longitude::Float64,
                      height::Int,
                      pointings::Array{Float64,1})

    τ_s = 1 / sampling_rate
    times = 0:τ_s:time_duration

    skydirs = Array{Float64}(undef, length(times), 2)
    skyψ = Array{Float64}(undef, length(times), 1)

    (skydirs, skyψ) = Sl.genpointings(pointings,
                                      times,
                                      t_start,
                                      t_stop;
                                      latitude_deg=latitude,
                                      longitude_deg=longitude,
                                      height_m=height) do time_s
                                          return (0.0, deg2rad(elevation_angle),
                                                  Sl.timetorotang(time_s,
                                                                  spin_velocity))
                                      end
    
    hits_map = Healpix.Map{Float64, Healpix.RingOrder}(NSIDE)
    for i in eachindex(skyψ)
        hits_map.pixels[Healpix.ang2pix(hits_map, π/2 - skydirs[i, 1],
                                        skydirs[i, 2])] += 1.0
    end
    hits_map
end

pars_arg = ARGS

NSIDE = 256

t_start = DateTime(2020, 04, 01, 0, 0, 0)
t_stop = DateTime(2022, 04, 01, 0, 0, 0)
time_duration_s = ((t_stop - t_start) / 1000).value 

tenerife_latitude_deg = Sl.TENERIFE_LATITUDE_DEG 
tenerife_longitude_deg = Sl.TENERIFE_LONGITUDE_DEG
tenerife_height_m = Sl.TENERIFE_HEIGHT_M 

elevation_angle_deg = parse(Float64, pars_arg[1])
spin_velocity_rpm = 1.0
sampling_rate_hz = 50.0

db = Sl.InstrumentDB()
keys_dbfp = collect(keys(db.focalplane))

for (j, f) in enumerate(keys_dbfp)
    if startswith(f, "W")
        delete!(db.focalplane, f)
        deleteat!(keys_dbfp, j)
    end
end

if size > length(keys_dbfp)
    N = Sl.split_into_n(length(keys_dbfp), length(keys_dbfp))
else
    N = Sl.split_into_n(length(keys_dbfp), size)
end

for r in rank
    z = r + 1
    if z > length(keys_dbfp)
    else
        subkeys = keys_dbfp[sum(N[1:z])-N[z]+1:sum(N[1:z])]
        for k in subkeys
            hits_map = get_hits_map(elevation_angle_deg,
                                    spin_velocity_rpm,
                                    sampling_rate_hz,
                                    time_duration_s,
                                    t_start,
                                    t_stop,
                                    NSIDE,
                                    tenerife_latitude_deg,
                                    tenerife_longitude_deg,
                                    tenerife_height_m,
                                    db.focalplane[k].orientation)
            
            save_to_path = join([pars_arg[2],
                                 db.focalplane[k].name,
                                 "_",
                                 Int(elevation_angle_deg),
                                 "°.fits"])
            
            Healpix.saveToFITS(hits_map, save_to_path, typechar = "D")
            
            f = FITS(save_to_path, "r+")
            
            write_key(f[2], "LATITUDE", tenerife_latitude_deg,
                      "Tenerife Latitude [Deg]")
            write_key(f[2], "LONGITUD", tenerife_longitude_deg,
                      "Tenerife Longitude [Deg]")
            write_key(f[2], "HEIGHT", tenerife_height_m, "Tenerife Height [m]")
            
            write_key(f[2], "DURATION", time_duration_s, "[sec]")
            write_key(
                f[2], "TSTART", Dates.format(t_start, "yyyy-mm-dd HH:MM:SS"),
                "Starting Date and time of the observation: yyyy-mm-dd HH:MM:SS")
            write_key(
                f[2], "TSTOP", Dates.format(t_stop, "yyyy-mm-dd HH:MM:SS"),
                "Ending Date and time of the observation: yyyy-mm-dd HH:MM:SS")
            
            write_key(f[2], "ANG_EL", elevation_angle_deg,
                      "Elevation angle [Deg]")
            write_key(f[2], "SPIN_VEL", spin_velocity_rpm, "Spin velocity [rpm]")
            write_key(f[2], "S_RATE", sampling_rate_hz, "Sampling rate [Hz]")
            
            write_key(f[2], "DETECTOR", db.focalplane[k].name, "Detector name")
            
            close(f)
        end
    end
end
    
MPI.Finalize()
    
