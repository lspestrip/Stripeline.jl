import Stripeline
const Sl = Stripeline

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
# The same as above, but with a 2× gain and an offset of 5 V

adc = Sl.ADC(offset_v = 5.0, gain_v_over_adu = 2.0)

@test Sl.adc_response(adc, 0.0) == -10
@test Sl.adc_response(adc, 1.0) == -8
@test Sl.adc_response(adc, -1.0) == -12
@test Sl.adc_response(adc, 2^21) == 2^19

@test Sl.adc_inv_response(adc, 0) ≈ 5.0
@test Sl.adc_inv_response(adc, 1) ≈ 5.5
@test Sl.adc_inv_response(adc, -1) ≈ 4.5
