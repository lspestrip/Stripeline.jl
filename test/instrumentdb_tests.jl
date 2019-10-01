import Stripeline
const Sl = Stripeline

db = Sl.InstrumentDB(joinpath(dirname(pathof(Stripeline)), "..", "instrumentdb"))

@test length(db.focalplane) ≥ 49
@test length(db.detectors) ≥ 55
@test "I0" in keys(db.focalplane)

# Now test that everything works without specifying the paths

defaultdb = Sl.InstrumentDB()

@test sort(collect(keys(defaultdb.focalplane))) == sort(collect(keys(db.focalplane)))
@test sort(collect(keys(defaultdb.detectors))) == sort(collect(keys(db.detectors)))

# Check the high-level API

@test Sl.detector(defaultdb, "I0") != nothing
@test Sl.detector(defaultdb, 2) != nothing

@test Sl.spectrum(defaultdb, 2) != nothing
@test Sl.spectrum(defaultdb, "I0") != nothing

@test Sl.tnoise(defaultdb, 2) != nothing
@test Sl.tnoise(defaultdb, "I0") != nothing

@test Sl.bandpass(defaultdb, 2) != nothing
@test Sl.bandpass(defaultdb, "I0") != nothing

x, y = Sl.bandshape(defaultdb, "I0")
@test length(x) > 0
@test length(x) == length(y)

f1_i, f1_q, f1_u = Sl.fknee_hz(defaultdb, 10)
f2_i, f2_q, f2_u = Sl.fknee_hz(defaultdb, 10, tsys_k = 1000.0)

@test f2_q < f1_q
@test f2_u < f1_u

# Just a few basic checks, to verify that the function does not crash
(sensitivity, num) = Sl.sensitivity_tant(defaultdb, 10.0)
@test sensitivity > 0.0
@test num == 49

(sensitivity, num) = Sl.sensitivity_tant(defaultdb, 10.0, modules = Set([0, 1]))
@test sensitivity > 0.0
@test num == 14

@test Sl.t_to_trj(2.72548, 43e9) ≈ 1.82263012068
@test Sl.trj_to_t(Sl.t_to_trj(2.72548, 43e9), 43e9) ≈ 2.72548
