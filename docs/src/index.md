```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Stripeline User's Manual

An implementation of a simulation/data analysis pipeline for the
LSPE/STRIP instrument.

To install Stripeline, start Julia and type the following command:
```julia
using Pkg
Pkg.add("https://github.com/lspestrip/Stripeline")
```

To run the test suite, type the following command:
```julia
using Pkg; Pkg.test("Stripeline")
```

In this manual, we will often assume that Stripeline has been imported
using the following commands:

```julia
import Stripeline
const Sl = Stripeline
```

In this way, we can call functions like [`genpointings`](@ref) using
the syntax `Sl.genpointings`, instead of the longer
`Stripeline.genpointings`.

## Documentation

The documentation was built using
[Documenter.jl](https://github.com/JuliaDocs).

```@example
println("Documentation built $(now()) with Julia $(VERSION).") # hide
```

## Index

```@index
```
