```@meta
DocTestSetup = quote
    using Stripeline
end
```
# Simulation tutorial

This aim of this tutorial is to describe how you can use the functions of Stripeline repository to perform a complete simulation of the LSPE/STRIP experiment. 

So far, the simulation includes:

- simulation of the telescope scanning the sky

- simulation of noise (white and 1/f)

- production of a tod 

- production of a map by using the destriping technique

Two examples will be presented: 

1. a simple case, in which we produce and analyze a very small TOD (1 day observation)

2. a more general and realistic case, suitable also for the production of large TODs.
In this case it will be necessary to use MPI functions.


## 1. Simple Case: small TOD

First of all, you should import the following packages:

```
import Healpix
imprt Random
import CorrNoise
import Stripeline
const Sl = Stripeline
using FITSIO

```
In this simple case we can avoid using MPI, but we need the following line: 
```
comm = missing

```
since some of the functions we will use (e.g. the `destripe` function) require the MPI communicator as input parameter. 



Then, you can start by setting the simulation parameters:

- the number of days of observation

- the number of polarimeters to simulate

- the sampling frequency (Hz)

- the length of 1/f baselines (s) as input for the destriper

- the NSIDE you prefer for your output map

- the temperature of the sky signal (K)

- the temperature of the atmosphere (K)

- the temperature of the telescope (K)

- the noise temperature of the polarimeters (K)

- the knee frequency of the polarimeters (Hz)

- the bandwidth of the polarimeters (Hz)

N.B. In the current version of the simulation we consider all polarimeters identical in properties.

```
#Simulation parameters

days_of_observation = 1
num_of_polarimeters = 1
baseline_length_s = 10  
fsamp_hz = 40
NSIDE = 256

tcmb_k = 3
tatm_k = 15
ttel_k = 3
tnoise_k = 35
fknee_hz = 0.01
β_hz = 7e9

```

At this point, you can compute some parameters we will need later: the total observation time (in s), the integration time, the total system temperature and the receiver sensitivity. 

```
total_time = days_of_observation * 24 * 3600
tsys_k = tnoise_k + tatm_k + ttel_k + tcmb_k
τ_s = 1 / fsamp_hz
σ_k = (tsys_k / sqrt(β_hz * τ_s))
```

Open now the input map, that is to say the sky you want to scan with your instrument. 
You can find two example input maps (with two different resolutions) here: https://github.com/lspestrip/Stripeline.jl/tree/master/test/testfiles.

Those maps have been produced with [PySM](https://github.com/bthorne93/PySM_public) and are specific for the STRIP case: they are 43 GHz maps of polarized emission only (cmb, synchrotron and dust).
If you want to produce your own input map, you can use this python [script](https://github.com/silviacaprioli/PySMforSTRIP/blob/master/PySMmap_production.py). 

```
inputmap = Healpix.readMapFromFITS(raw"PySM_inputmap_nside256.fits", 2 , Float64)
inputmap_resol = inputmap.resolution

resol = Healpix.Resolution(NSIDE) #desired resolution for output map
num_of_pixels = resol.numOfPixels
```
You can now scan the input map according to your scanning strategy, thus producing the pure signal "sky TOD".

In this case we use nominal scanning strategy for STRIP (Tenerife latitude, 20 degrees of elevation angle, 1 rpm)

```
#Generate sky tod

times = 0:τ_s: (total_time-τ_s)
(dirs, ψ) = Sl.genpointings([0, 0, 1], times; latitude_deg=28.29) do time_s
    return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, 1))
end

pix_idx_inputmap = Healpix.ang2pixRing.(Ref(inputmap_resol), dirs[:, 1], dirs[:, 2])
pix_idx = Healpix.ang2pixRing.(Ref(resol), dirs[:, 1], dirs[:, 2])

sky_tod = inputmap.pixels[pix_idx_inputmap]
```

Now you should add noise to your sky tod. 

We simulate both white noise and 1/f noise, in accordance with the noise properties of the polarimeters specified at the beginning of the script.

To do that, we use the functions of the module [CorrNoise.jl](https://github.com/ziotom78/CorrNoise.jl).
```
#Generate Noise

seed = rand(1:1000)
rng = CorrNoise.OofRNG(CorrNoise.GaussRNG(Random.MersenneTwister(seed)), -1, 1.15e-5, fknee_hz, fsamp_hz);
noise_tod = [CorrNoise.randoof(rng) * σ_k for i in 1:(fsamp_hz * total_time)]

```
and, finally, you can get the final, realistic TOD just by doing:
```
tod = sky_tod + noise_tod
```
Once simulated your data, you can now perform data analysis.

You can call the destriper and clean the map from 1/f noise.

(N.B. the destriper needs in input an array containing the lengths of all 1/f baselines.
For the sake of semplicity, we consider baselines of equal length).

```
#Run the destriper

baseline_len = repeat([baseline_length_s*fsamp_hz],Int64(total_time/baseline_length_s))

(destr_map, a) = Sl.destripe(pix_idx, tod, num_of_pixels, baseline_len, comm)
```
The output of the destriper are:

- `destr_map` : the destriped map, cleaned from 1/f noise.
- `a` : the baselines array.


If you wish, you can finally save the destriped map in a .fits file: 
```
#save file 

mapfile = Healpix.Map{Float64,Healpix.RingOrder}(NSIDE)
mapfile.pixels = destr_map
Healpix.saveToFITS(mapfile, "!results/destriped_map.fits", typechar = "D")
```
