```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Instrument database

`InstrumentDB` takes advantage of the structures [`Detector`](@ref) and
[`Horn`](@ref) to retrieve information about feed horns and detectors from
a YAML file. There are a set of YAML files containing the default configuration
for the STRIP instrument in the repository.

The following example initializes an object of type [`InstrumentDB`](@ref) with
the values referred to the standard STRIP instrument:

```@repl instrumentdbexample
using Stripeline; # hide
db = InstrumentDB()
```

As `db` is a `struct`, its field can be accessed with the usual dot notation. The
two fields in `db` are `focalplane` and `detectors`. They are both dictionaries,
associating horn names to `Horn` objects and detectors IDs to `Detector` objects,
respectively.


```@docs
InstrumentDB
```

```@docs
Horn
```

```@docs
Detector
BandshapeInfo
SpectrumInfo
NoiseTemperatureInfo
```