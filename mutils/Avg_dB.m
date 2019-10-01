function [ avg ] = Avg_dB(dB, dim, coef)
%AVG_DB This function calculates true average of variable in dB.
%   The function calculates average of absolute variable given in dB along
%   dimension dim. Coef is defined in original equation 
%   dB = coef * log10(a). It is 10 for power and 20 for voltage.
%   The returned value is again in dB.

% Authors: Jan Kral <kral.j@lit.cz>
% Date: 12.1.2017

avg = coef*log10(sum(10.^(dB/coef), dim)./size(dB,dim));

end

