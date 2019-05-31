# -*- encoding: utf-8 -*-

export ADC, optimized_adc, adc_response, adc_inv_response, adc_filter

@doc raw"""
A structure representing the configuration of an Analogue-to-Digital Converter
(ADC). To simplify things, we assume that the ADC is fed with some *temperature*,
instead of a voltage. In this case we can avoid dealing with the calibration from
Volt to Kelvin, which is seldom useful in simulations.

The equation used to model the ADC is the following:

    x(V) = round(gain_k_over_adu * (V - offset_k)) + zero_point_adu

and the output is clipped within a user-defined range. The following fields are
available:

- `offset_k`
- `gain_k_over_adu`
- `zero_point_adu`
- `min_output_adu`: output values below this will be clipped
- `max_output_adu`: output values above this will be clipped
- `non_linearities_x_adu`: see below
- `non_linearities_y_adu`: see below

Non linearities must be specified using the two arrays `non_linearities_x_adu`
and `non_linearities_y_adu`. For each input voltage fed to the ADC, a non
linearity is applied to the output according to the following algorithm:

- The function computes the ideal output (i.e., without non linearities);
- It checks if the value of the ideal output is found in the array
  `non_linearities_x_adu`;
- If the value is found, add the corresponding value in the array
  `non_linearities_y_adu` (i.e., the element with the same index as the element
  in `non_linearities_x_adu`) to the ideal output;
- If the value is not found, but it falls within two consecutive values,
  use a linear interpolation;
- If the value is smaller than the first element or larger than the last element
  in `non_linearities_x_adu`, do not apply any correction to the output.

It is **fundamental** that `non_linearities_x_adu` is sorted in ascending
order, and that the order of the in elements `non_linearities_y_adu` matches
the order in `non_linearities_x_adu`.

"""
mutable struct ADC
    offset_k::Float64
    gain_k_over_adu::Float64
    zero_point_adu::Float64
    min_output_adu::Int
    max_output_adu::Int
    non_linearities_x_adu::Array{Float64,1}
    non_linearities_y_adu::Array{Float64,1}

    ADC(;
        offset_k = 0.0,
        gain_k_over_adu = 1.0,
        zero_point_adu = 0.0,
        min_output_adu = -2^19,
        max_output_adu = 2^19,
        non_linearities_x_adu = Float64[],
        non_linearities_y_adu = Float64[]) = new(offset_k, gain_k_over_adu, zero_point_adu, min_output_adu, max_output_adu, non_linearities_x_adu, non_linearities_y_adu)
end

@doc raw"""

    optimized_adc(; min_input_k = 0.0, max_input_k = 100.0, dynamic_range = 0.35, nbits = 20, non_linearities_x_adu = Float64[], non_linearities_y_adu = Float64[])

Return an object of type `ADC` that is optimized to measure temperatures between
`min_input_k` and `max_input_k`. The ADC is configured to return *signed*
numbers in the range -2^nbits…2^nbits.

The function accepts the following keywords:

- `min_input_k`: temperature that should trigger the lowest output within the
  dynamic range
- `max_input_k`: temperature that should trigger the lowest output within the
  dynamic range
- `dynamic_range`: pure number in the interval 0…1 that specifies the dynamic
  range of the output
- `nbits`: number of bits used by the ADC
- `non_linearities_x_adu`: same parameter used in the constructor for `ADC`
- `non_linearities_y_adu`: same parameter used in the constructor for `ADC`

"""
function optimized_adc(; min_input_k = 0.0, max_input_k = 100.0, dynamic_range = 0.35, nbits = 20, non_linearities_x_adu = Float64[], non_linearities_y_adu = Float64[])
    min_output_adu = -2^(nbits - 1)
    max_output_adu = 2^(nbits - 1)
    ADC(offset_k = -min_input_k,
        gain_k_over_adu = (max_input_k - min_input_k) / (2^nbits),
        zero_point_adu = min_output_adu,
        min_output_adu = min_output_adu,
        max_output_adu = max_output_adu,
        non_linearities_x_adu = non_linearities_x_adu,
        non_linearities_y_adu = non_linearities_y_adu,
    )
end

@doc raw"""
    adc_response(adc::ADC, input_k; include_non_linearities = true)

Simulate the response of an ADC on some voltage `input_k`. If
`include_non_linearities` is `true`, non-linearities specified in `adc` will be
considered; otherwise, the ADC will be assumed to be ideal.

See also [`adc_inv_response`](@ref) for the (pseudo)inverse function.
"""
function adc_response(adc::ADC, input_k; include_non_linearities = true)
    output = adc.gain_k_over_adu * (input_k - adc.offset_k) + adc.zero_point_adu
    
    # Let's use some shorthands here
    x = adc.non_linearities_x_adu
    y = adc.non_linearities_y_adu

    if include_non_linearities && length(x) > 0 && (x[1] ≤ output ≤ x[end])
        idx = searchsortedlast(x, output)
        if idx > 0
            δ = if idx == length(x)
                y[idx]
            else
                # Interpolate between idx and idx + 1
                y[idx] + (output - x[idx]) * (y[idx] - y[idx + 1]) / (x[idx] - x[idx + 1])
            end
            output += δ
        end
    end

    clamp(round(Int, output),
        adc.min_output_adu, 
        adc.max_output_adu,
    )
end

@doc raw"""
    adc_inv_response(adc::ADC, input_adu)
    
Apply the inverse transformation of an ADC to get some voltage from a digital
measurement in ADU.This function assumes that the ADC is ideal, i.e., that it
does not have non-linearities.

See also [`adc_response`](@ref).
"""
function adc_inv_response(adc::ADC, input_adu)
    (input_adu - adc.zero_point_adu) / adc.gain_k_over_adu + adc.offset_k
end

@doc raw"""
    adc_filter(adc::ADC, input_k; include_non_linearities = true)

Simulate the measurement of the temperature `input_k` through the ADC `adc`. The
result is still a temperature, after it has been converted into ADUs by the ADC
and then converted back to a temperature again. If non-linearities are specified
in `adc`, they will be applied only to the K→ADU transformation.
"""
function adc_filter(adc::ADC, input_k; include_non_linearities = true)
    adu = adc_response(adc, input_k, include_non_linearities = include_non_linearities)
    adc_inv_response(adc, adu)
end