function [Idiff_adc,I_res,res_bit] = ADC(Idiff,min_I,max_I,n_bit_adc,max_adc_bit)
signIdiff = sign(Idiff);
Idiff_abs = abs(Idiff);
if ~exist('max_adc_bit', 'var')
    max_adc_bit = 12;
end
res_bit = 0;
if n_bit_adc> max_adc_bit
    res_bit = n_bit_adc-max_adc_bit;
    n_bit_adc = max_adc_bit;
end
Imax = max(Idiff_abs,[],"all");

I_res = (max_I-min_I)/(2^n_bit_adc-1);
% assert(Imax<=max_I+5*I_res)
Idiff_abs(Idiff_abs>max_I) = max_I;
Idiff_abs(Idiff_abs<min_I) = min_I;
% Idiff_adc = signIdiff.*(round(Idiff_abs/I_res)*I_res);
Idiff_adc_tmp = signIdiff.*(round(Idiff_abs/I_res));
Idiff_adc = Idiff_adc_tmp*2.^(res_bit);
end