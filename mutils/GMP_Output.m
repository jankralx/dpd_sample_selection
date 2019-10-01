function y = GMP_Output(x, coefs, Ka, La, Kb, Lb, Mb, Kc, Lc, Mc, type)

% GMP_Output function returns output of GMP model create from input signal
%   x. Generalized Memory Polynomial (GMP) model has nonlinearity order Ka,
%   Kb, Kc, memory length La, Lb, Lc, and number of cross-terms Mb, Mc.
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the DDR2 matrix U
%   
%   coefs - coefficients of GMP model
%
%   type == 'num' (default):
%
%   Ka, Kb, Kc  - maximum DPD nonlinearity order
%
%   La, Lb, Lc  - maximum DPD memory length
%
%   Mb, Mc - number of cross-terms
%
%   type == 'all':
%   
%   Ka, Kb, Kc  - vector of DPD nonlinearity orders evaluated
%
%   La, Lb, Lc  - vector of DPD memory length evaluated
%
%   Mb, Mc - vector of cross-terms evaluated
%
%   Returns:
%   ========
%
%   y  - output of GMP model with coefficients coefs (otype = 'modelout').


% Authors: Jan Kral <kral.j@lit.cz>
% Date: 10.2.2019

if nargin < 11
    if length(Ka) > 1 || length(La) > 1 || length(Kb) > 1 || ...
            length(Lb) > 1 || length(Mb) > 1 || length(Kc) > 1 || ...
            length(Lc) > 1 || length(Mc) > 1
        type = 'all';
    else
        type = 'num';
    end
end

y = GMP_Matrix(x, Ka, La, Kb, Lb, Mb, Kc, Lc, Mc, type, 'modelout', coefs);
