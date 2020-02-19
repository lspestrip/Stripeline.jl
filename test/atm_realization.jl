# To do...
import Healpix
import Random
import CorrNoise
import Stripeline
const Sl = Stripeline
using FITSIO
using Plots

N = 252
cell_size = 0.2 # km/pix
corr_length = 0.5 #km

rng, random_realization, mag_k = Sl.initialize_kolmogorov(N, 123, cell_size, corr_length)

atm = Sl.atm_observe(random_realization, mag_k ./ 0.0001, 0.0, 70.0)

grad = ColorGradient([:blue, :white])
heatmap(atm[:,:,10], c=grad)

for i in 1:10
    ran_real = Sl.atm_update(random_realization, rng, N, 1.0, 1.0, cell_size, 50.0)
    atm = Sl.atm_observe(ran_real, mag_k, 0.0, 70.0)
    heatmap(atm[:,:,100], c=grad)
    savefig(string(i))
end
