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

A number of high-level functions ease the access of the fields in a
[`InstrumentDB`](@ref) object:

- [`detector`](@ref) returns a [`Detector`](@ref) structure, containing the
  details of a polarimeter;
- [`bandpass`](@ref) returns a [`BandshapeInfo`](@ref) structure,
  containing the shape of the bandpass of a detector;
- [`spectrum`](@ref) returns a [`SpectrumInfo`](@ref)
- [`fknee`](@ref) returns the knee frequency of the 1/f noise for the I, Q, and
  U signals, adapted to the brightness temperature of the load being observed by
  the detector;
- [`tnoise`](@ref) returns the noise temperature for the I, Q, and U components.

The structure `Detector` uses three structures to organize its data in a
hierarchical way:

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
bandshape
SpectrumInfo
NoiseTemperatureInfo
```

## High-level access functions

```@docs
detector
bandpass
spectrum
fknee_hz
tnoise
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
