db = InstrumentDB(joinpath(Pkg.dir("Stripeline"), "instrumentdb"))

@test length(db.focalplane) ≥ 49
@test length(db.detectors) ≥ 55