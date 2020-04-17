```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Basic functions

## Measure unit conversion

It is often useful to convert measurements between thermodynamic
temperatures and Rayleigh-Jeans temperatures. The following functions
implement this kind of conversion. Note that there are two families of
functions:

1. Functions that convert *absolute* measurements: [`t_to_trj`](@ref),
   [`trj_to_t`](@ref);
2. Functions that convert *sensitivities* (i.e., small fluctuations
   around an absolute value): [`deltat_to_deltatrj`](@ref),
   [`deltatrj_to_deltat`](@ref)

The function [`sensitivity_tant`](@ref) computes the overall
sensitivity of a set of polarimeters, using information from [The
Strip Instrument Database](@ref).

```@docs
sensitivity_tant
t_to_trj
trj_to_t
deltat_to_deltatrj
deltatrj_to_deltat
```
