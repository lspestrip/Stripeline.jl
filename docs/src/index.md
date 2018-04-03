```@meta
DocTestSetup = quote
    using Stripeline2
end
```

# Stripeline User's Manual

An implementation of a simulation/data analysis pipeline for the LSPE/STRIP instrument.

To install Stripeline2, start Julia and type the following command:
```julia
Pkg.clone("https://github.com/lspestrip/stripeline2")
```

To run the test suite, type the following command:
```julia
Pkg.test("stripeline2")
```

In this manual, we will often assume that Stripeline2 has been imported using the following commands:
```julia
import Stripeline2
const Sl = Stripeline2
```
In this way, we can call functions like [`genpointings`](@ref) using the syntax `Sl.genpointings`, instead of the longer `Stripeline2.genpointings`.

## Documentation

The documentation was built using [Documenter.jl](https://github.com/JuliaDocs).

```@example
println("Documentation built $(now()) with Julia $(VERSION).") # hide
```

## Index

```@index
```
