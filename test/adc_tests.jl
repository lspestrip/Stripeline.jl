import Stripeline

##########################################################################
# Simple ADC, with gain 1, no offset and 20 bit of dynamic range

adc = Sl.ADC()

@test Sl.adc_response(adc, 0.0) == 0
@test Sl.adc_response(adc, 1.0) == 1
@test Sl.adc_response(adc, -1.0) == -1
@test Sl.adc_response(adc, 2^21) == 2^19

@test Sl.adc_inv_response(adc, 0) ≈ 0.0
@test Sl.adc_inv_response(adc, 1) ≈ 1.0
@test Sl.adc_inv_response(adc, -1) ≈ -1.0

##########################################################################
# The same as above, but with a 2× gain and an offset of 5 K

adc = Sl.ADC(offset_k = 5.0, gain_k_over_adu = 2.0)

@test Sl.adc_response(adc, 0.0) == -10
@test Sl.adc_response(adc, 1.0) == -8
@test Sl.adc_response(adc, -1.0) == -12
@test Sl.adc_response(adc, 2^21) == 2^19

@test Sl.adc_inv_response(adc, 0) ≈ 5.0
@test Sl.adc_inv_response(adc, 1) ≈ 5.5
@test Sl.adc_inv_response(adc, -1) ≈ 4.5

##########################################################################
# Non-linearity checks

adc = Sl.ADC(
    min_output_adu = -10,
    max_output_adu = 10,
    non_linearities_x_adu = [-10, 0, 10],
    non_linearities_y_adu = [0, 10, 0],
)

@test Sl.adc_response.(Ref(adc), [-10.0, -5.0, 0.0, 5.0, 10.0]) == [-10, 0, 10, 10, 10]

##########################################################################
# Additional checks

adc = Sl.ADC()

@test Sl.adc_filter(adc, 1.1) ≈ 1.0