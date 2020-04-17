```@meta
DocTestSetup = quote
    using Stripeline
end
```

# The Strip Instrument Database

## Introduction

Every simulator of a scientific experiment needs an «instrument
database», i.e., a way to access the specification of the instrument
to be simulated. The kind of information stored in this database is
the following:

-   How many detectors are present in the instrument, and what are
    their characteristics in terms of noise and placement on the focal
    plane;
-   What is the response of the optical elements, i.e., the so-called
    *beam function*;
-   Which frequencies can be measured by each detector (the so-called
    *bandpass*);
-   Ect.

Stripeline implements an instrument database using [YAML
files](https://en.wikipedia.org/wiki/YAML) and two data structures,
[`Detector`](@ref) and [`Horn`](@ref), plus some functions to easily
access the information in them.

The following example initializes an object of type
[`InstrumentDB`](@ref) with the values referred to the standard STRIP
instrument:

```@repl instrumentdbexample
using Stripeline; # hide
db = InstrumentDB()
```

This command loads the YAML files provided in the Stripeline
repository and initializes the `db` object. As `db` is a `struct`, its
field can be accessed with the usual dot notation. The two fields in
`db` are `focalplane` and `detectors`. They are both dictionaries,
associating horn names to `Horn` objects and detectors IDs to
`Detector` objects, respectively:

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
- [`fknee_hz`](@ref) returns the knee frequency of the 1/f noise for the I, Q, and
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

Some of them have the ability to be plotted using
[Plots.jl](https://github.com/JuliaPlots/Plots.jl):

```@repl instrumentdbexample
using Plots
plot(db.detectors[2].bandshape)
savefig("instrumentdb-bandshape.svg"); nothing # hide
```

![](instrumentdb-bandshape.svg)

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
