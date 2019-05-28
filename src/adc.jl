# -*- encoding: utf-8 -*-

export ADC, adc_response, adc_inv_response

@doc raw"""
A structure representing the configuration of an Analogue-to-Digital Converter
(ADC). The equation used to model the ADC is the following:

    x(V) = round(gain_v_over_adu * (V - offset_v)) + zero_point_adu

and the output is clipped using a predefined number of bits. The following
fields are available:

- `offset_v`
- `gain_v_over_adu`
- `zero_point_adu`
- `min_output_adu`
- `max_output_adu`
- `non_linearities_x_adu`
- `non_linearities_y_adu`
"""
struct ADC
    offset_v::Float64
    gain_v_over_adu::Float64
    zero_point_adu::Float64
    min_output_adu::Int
    max_output_adu::Int
    non_linearities_x_adu::Array{Float64,1}
    non_linearities_y_adu::Array{Float64,1}

    ADC(;
        offset_v = 0.0,
        gain_v_over_adu = 1.0,
        zero_point_adu = 0.0,
        min_output_adu = -2^19,
        max_output_adu = 2^19,
        non_linearities_x_adu = Float64[],
        non_linearities_y_adu = Float64[]) = new(offset_v, gain_v_over_adu, zero_point_adu, min_output_adu, max_output_adu, non_linearities_x_adu, non_linearities_y_adu)
end


@doc raw"""
    adc_response(adc::ADC, input_v; include_non_linearities = true)

Simulate the response of an ADC on some voltage `input_v`. If
`include_non_linearities` is `true`, non-linearities specified in `adc` will be
considered; otherwise, the ADC will be assumed to be ideal.

See also [`adc_inv_response`](@ref) for the (pseudo)inverse function.
"""
function adc_response(adc::ADC, input_v; include_non_linearities = true)
    output = adc.gain_v_over_adu * (input_v - adc.offset_v) + adc.zero_point_adu
    
    if include_non_linearities
        # Let's use some shorthands here
        x = adc.non_linearities_x_adu
        y = adc.non_linearities_y_adu

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
    (input_adu - adc.zero_point_adu) / adc.gain_v_over_adu + adc.offset_v
end
