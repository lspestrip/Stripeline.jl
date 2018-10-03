using Test
import Stripeline

eps = 5e-5 

t_start = "2018-01-01 00:00:00"
t_stop = "2018-01-01 00:00:00"

crab_az_stellarium_rad = 3.168070267969909
crab_alt_stellarium_rad = 1.4612458881566615

dirs = [π/2-crab_alt_stellarium_rad 2π-crab_az_stellarium_rad]

crab_ra_astropy_rad = 1.4596726619436968   
crab_dec_astropy_rad = 0.3842255081802917 
crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

skydirs = genskypointings(t_start,
                          t_stop,
                          dirs;
                          latitude_deg=TENERIFE_LATITUDE_DEG,
                          longitude_deg=TENERIFE_LONGITUDE_DEG,
                          height_m=TENERIFE_HEIGHT_M)
crab_position_skydirs = sqrt(skydirs[1, 1]^2 + skydirs[1, 2]^2)

@test skydirs[1, 1] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_dec_astropy_rad atol = eps
@test crab_position_skydirs ≈ crab_position atol = eps
