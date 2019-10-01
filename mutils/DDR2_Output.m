function y = DDR2_Output(x, K, M, coef)

% DDR2_Output function returns signal y created by
%   transformation of signal x by DDR2 model
%   with maximum order K and memory length M. Model is
%   characterised by coeficients coef.
%
%   Inputs:
%   =======
%   x  - signal to be transformed by DDR2 model
%
%   K  - nonlinearity order
%
%   M  - memory length
%
%   coef - coefficients of the MP model
%
%   Returns:
%   ========
%
%   y - output of DDR2 model - transformed signal x

% Authors: Jan Kral <kral.j@lit.cz>, Tomas Gotthans
% Date: 19.8.2018

y = DDR2_Matrix(x, K, M, 'modelout', coef);