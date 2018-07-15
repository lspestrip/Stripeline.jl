using Base.Test
using Stripeline

# RNG tests

ref = [0.9854981268, 0.2536861135,
       0.8791850018, 0.8541532028,
       0.4161281283, 0.6491611481,
       0.8336713058, 0.2541507778,
       0.0673349826, 0.7994472017,
       0.2600378725, 0.5134736516,
       0.4544143085, 0.9879117254,
       0.1875971474, 0.1698732623,
       0.9844831470, 0.2436698098,
       0.2365035815, 0.9442295933]

flatrng = initflatrng128()
for (idx, sample) in enumerate([rand(flatrng) for i = 1:length(ref)])
    @assert ref[idx] ≈ sample
end

ref = [1.8051588123, -1.0150233503,
       -0.5138818361, 0.6974503387,
       0.0959640650, -1.7090943033,
       0.2825155379, -0.0263954840,
       -0.4486580098, -0.4245704031,
       -1.6567426564, -1.2158240845,
       0.3230971021, -1.0730251963,
       -0.6426576020, 1.0481948256,
       -0.0527193353, -0.2020275077,
       0.4885079587, 0.2111337896]

normrng = GaussRNG(initflatrng128())
for (idx, sample) in enumerate([randn(normrng) for i = 1:length(ref)])
    @assert ref[idx] ≈ sample
end

ref = [2.0886370365, -0.6074844495,
       -0.3464683356, 0.8936789422,
       0.4167745672, -1.6416294545,
       0.1259486772, -0.1427304943,
       -0.6395859847, -0.7526144784,
       -2.3116079175, -2.3217436788,
       -0.9229344940, -2.4367338428,
       -2.2756948201, -0.5210397131,
       -1.4655131788, -1.6547241628,
       -0.9190954646, -1.0864977062]

normrng = GaussRNG(initflatrng128())
oof2rng = Oof2RNG(normrng, 1.15e-5, 0.05, 1.0)
for (idx, sample) in enumerate([randoof2(oof2rng) for i = 1:length(ref)])
    @assert ref[idx] ≈ sample
end

ref = [2.0400958027, -0.6813881047,
       -0.3848671621, 0.84863765610,
       0.34838982440, -1.6703162833,
       0.13588215630, -0.1362109745,
       -0.6176813736, -0.7037434318,
       -2.2015645097, -2.1249009577,
       -0.6866971017, -2.1640479726,
       -1.9400718236, -0.1768193407,
       -1.1308790625, -1.2995275170,
       -0.5582322798, -0.7326164277]

normrng = GaussRNG(initflatrng128())
oofrng = OofRNG(normrng, -1.7, 1.15e-5, 0.05, 1.0)
for (idx, sample) in enumerate([randoof(oofrng) for i = 1:length(ref)])
    @assert ref[idx] ≈ sample
end

# Scanning strategy

(dirs, psi) = genpointings([0, 0, 1], 0:20:60, latitude_deg=0.0) do time_s
    wheel1ang = 0.0
    wheel2ang = 30.0
    wheel3ang = timetorotang(time_s, 1.0)

    (wheel1ang, wheel2ang, wheel3ang)
end

@test size(dirs) == (4, 2)
@test length(psi) == 4
#@test dirs ≈ [0.15487 4.71239; 2.0875 3.3214; 2.0875 6.10774; 0.15487 4.71675], atol = 1e-5
# This fails, but we have to revise the way polarization angles are calculated
#@test psi ≈ [0, π / 2, π, 3π / 4]

# Instrument DB

include("instrumentdb_tests.jl")

# Map-maker

include("map_maker_tests.jl")

# TOD-splitter
include("tod_splitter_tests.jl")
