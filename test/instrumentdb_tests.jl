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

# Just a few basic checks, to verify that the function does not crash
(sensitivity, num) = Sl.sensitivity_tant(defaultdb, 10.0)
@test sensitivity > 0.0
@test num == 49

(sensitivity, num) = Sl.sensitivity_tant(defaultdb, 10.0, modules = Set([0, 1]))
@test sensitivity > 0.0
@test num == 14

@test Sl.t_to_trj(2.72548, 43e9) ≈ 1.82263012068
@test Sl.trj_to_t(Sl.t_to_trj(2.72548, 43e9), 43e9) ≈ 2.72548