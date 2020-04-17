```@meta
DocTestSetup = quote
    using Stripeline
end
```

# Simulating data acquisition

An essential part of Strip polarimeters is the set of four
Analogue-to-Digital Converters (ADC) that measure input voltages as
20-bit digital numbers. This process is called *quantization*, and it
is a non-invertible operation that causes loss of information. Ideal
ADCs perform a linear operation (modulo a rounding operation), but
real-world components are never perfectly linear.

Because of the fact that CMB experiments like Strip measure brightness
temperatures, Stripeline models ADCs as devices that convert
temperatures into ADUs, neglecting the fact that Strip polarimeters
convert incoming fluxes into voltages.

Stripeline offers a few functions to simulate the behaviour of an ADC.
The simulation of ADC behaviour is useful to estimate the impact of
quantization and non linearities. The following schema show how things
work:

![](assets/adc_functions.svg)


Function [`adc_response`](@ref) takes a temperature as input, and it
produces the output that would be emitted by an ADC. The function
[`adc_inv_response`](@ref) performs the reverse transformation: it
converts a digital number back to a temperature. Function
[`adc_filter`](@ref) combines the two functions: it takes a
temperature as input, and it returns the temperature that has been
measured by the ADC, including the effect of quantization and non
linearities.

An important difference between `adc_response` and `adc_inv_response`
is the fact that `adc_response` considers non linearities, while
`adc_inv_response` does not.

```@docs
ADC
optimized_adc
adc_response
adc_inv_response
adc_filter
```
