import Stripeline

db = InstrumentDB(joinpath(dirname(pathof(Stripeline)), ".." ,"instrumentdb"))

@test length(db.focalplane) ≥ 49
@test length(db.detectors) ≥ 55
@test "I0" in keys(db.focalplane)

# Now test that everything works without specifying the paths

defaultdb = InstrumentDB()

@test sort(collect(keys(defaultdb.focalplane))) == sort(collect(keys(db.focalplane)))
@test sort(collect(keys(defaultdb.detectors))) == sort(collect(keys(db.detectors)))