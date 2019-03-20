```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Instrument database

`InstrumentDB` takes advantage of the structures [`Detector`](@ref)
and [`Horn`](@ref) to retrieve information about feed horns and
detectors from a YAML file. There are a set of YAML files containing
the default configuration for the STRIP instrument in the repository.

## Quick introduction

The following example initializes an object of type
[`InstrumentDB`](@ref) with the values referred to the standard STRIP
instrument:

```@repl instrumentdbexample
using Stripeline; # hide
db = InstrumentDB()
```

As `db` is a `struct`, its field can be accessed with the usual dot
notation. The two fields in `db` are `focalplane` and
`detectors`. They are both dictionaries, associating horn names to
`Horn` objects and detectors IDs to `Detector` objects, respectively:

```@repl instrumentdbexample
db.focalplane["I0"]
db.detectors[2]
```

The structure `Detector` is complex, as it is built over three other
structures:

- [`BandshapeInfo`](@ref)
- [`SpectrumInfo`](@ref)
- [`NoiseTemperatureInfo`](@ref)

All these structures know how to show themselves on the REPL:

```@repl instrumentdbexample
db.detectors[2].bandshape
db.detectors[2].spectrum
db.detectors[2].tnoise
```

For more information about the fields in the structures listed above,
as well as their meaning, keep reading.

## Structures

```@docs
InstrumentDB
Horn
Detector
BandshapeInfo
SpectrumInfo
NoiseTemperatureInfo
```

## Loading custom databases

It is not needed to load the default instrument database, as
Stripeline provides a number of additional functions to build mock
databases from dictionaries.

```@docs
defaultdbfolder
parsefpdict
parsedetdict
```
