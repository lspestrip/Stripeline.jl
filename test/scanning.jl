using Test
import Stripeline

eps = 1e-4 

t_start = "2018-01-01 00:00:00"
t_stop = "2018-01-01 00:00:00"

crab_az_stellarium_rad = 3.0925886563306033
crab_alt_stellarium_rad = 1.4612066182484915

dirs = [π/2-crab_alt_stellarium_rad 2π-crab_az_stellarium_rad]

crab_ra_rad = 1.4596726619436968 #83.633083° - astropy  
crab_dec_rad = 0.3842255081802917 #22.0145° - astropy

# crab_ra_rad = 1.4650141024930854 #83.939125 - stellarium
# crab_dec_rad = 0.3844310691810821 #22.02627777777778 - stellarium

skydirs = genskypointings(t_start,
                          t_stop,
                          dirs;
                          latitude_deg=TENERIFE_LATITUDE_DEG,
                          longitude_deg=TENERIFE_LONGITUDE_DEG,
                          height_m=TENERIFE_HEIGHT_M)

@test skydirs[1, 1] ≈ crab_ra_rad atol = eps
@test skydirs[1, 2] ≈ crab_dec_rad atol = eps
