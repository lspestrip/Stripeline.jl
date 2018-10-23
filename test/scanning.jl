using Test
import Stripeline

eps = 1e-4 


t_start = "2018-01-01 00:00:00"
t_stop = "2018-01-01 00:00:00"

crab_az_stellarium_rad = 3.168070267969909
crab_alt_stellarium_rad = 1.4612458881566615

dirs = [π/2-crab_alt_stellarium_rad 2π-crab_az_stellarium_rad]

crab_ra_astropy_rad = 1.4596726619436968   
crab_dec_astropy_rad = 0.3842255081802917 
crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

skydirs = genastropyskypointings(t_start,
                                 t_stop,
                                 dirs;
                                 latitude_deg=TENERIFE_LATITUDE_DEG,
                                 longitude_deg=TENERIFE_LONGITUDE_DEG,
                                 height_m=TENERIFE_HEIGHT_M)
crab_position_skydirs = sqrt(skydirs[1, 1]^2 + skydirs[1, 2]^2)

@test skydirs[1, 1] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_dec_astropy_rad atol = eps
@test crab_position_skydirs ≈ crab_position atol = eps

# skydirs = genskypointings(t_start,
#                           t_stop,
#                           dirs;
#                           latitude_deg=TENERIFE_LATITUDE_DEG,
#                           longitude_deg=TENERIFE_LONGITUDE_DEG,
#                           height_m=TENERIFE_HEIGHT_M)
# crab_position_skydirs = sqrt(skydirs[1, 1]^2 + skydirs[1, 2]^2)

# @test skydirs[1, 1] ≈ crab_ra_astropy_rad atol = eps
# @test skydirs[1, 2] ≈ crab_dec_astropy_rad atol = eps
# @test crab_position_skydirs ≈ crab_position atol = eps


t_start = "2019-02-01 02:00:00"
t_stop = "2019-02-03 02:00:00"

dirs = [π/2-0.6174703397894518 2π-4.852881197778698 ;
        π/2-0.6024799007695448 2π-4.859519751514131 ;
        π/2-0.5875049757874335 2π-4.866139397516001]

crab_ra_astropy_rad = 1.4596726619436968   
crab_dec_astropy_rad = 0.3842255081802917 
crab_position = sqrt(crab_ra_astropy_rad^2 + crab_dec_astropy_rad^2)

skydirs = genastropyskypointings(t_start,
                                 t_stop,
                                 dirs;
                                 latitude_deg=TENERIFE_LATITUDE_DEG,
                                 longitude_deg=TENERIFE_LONGITUDE_DEG,
                                 height_m=TENERIFE_HEIGHT_M)
crab_position_skydirs = sqrt.(skydirs[:, 1].^2 + skydirs[:, 2].^2)

@test skydirs[1, 1] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[1, 2] ≈ crab_dec_astropy_rad atol = eps
@test skydirs[2, 1] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[2, 2] ≈ crab_dec_astropy_rad atol = eps
@test skydirs[3, 1] ≈ crab_ra_astropy_rad atol = eps
@test skydirs[3, 2] ≈ crab_dec_astropy_rad atol = eps
@test crab_position_skydirs[1] ≈ crab_position atol = eps
@test crab_position_skydirs[2] ≈ crab_position atol = eps
@test crab_position_skydirs[3] ≈ crab_position atol = eps
