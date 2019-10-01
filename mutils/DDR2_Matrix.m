function U = DDR2_Matrix(x, K, M, otype, coefs)

% DDR2_matrix function returns matrix U created from the given signal x.
%   Matrix is created for Dynamic Deviation Reduction (DDR2) model
%   of nonlinearity order K and memory length M.
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the DDR2 matrix U
%
%   K  - DPD nonlinearity order
%
%   M  - DPD memory length
%
%   [otype] - 'matrix' : function returns matrix created from signal x
%                        (default)
%           - 'modelout' : function returns output signal of GMP model with
%                          coefficients coefs
%
%   [coefs] - coefficients used for signal output calculation
%
%   Returns:
%   ========
%
%   U  - DDR2 matrix created from signal x (otype = 'matrix') or output
%        signal of GMP model with coefficients coefs (otype = 'modelout').

% Authors: Jan Kral <kral.j@lit.cz>, Tomas Gotthans
% Date: 21.1.2018
% Update: 22.3.2019 - speed optimization


% the DDR2 is defined in the article
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5773325 
% Equation 3

if nargin < 4
    otype = 'matrix';
end

% convert otype to numeric for faster comparison
switch otype
    case 'matrix'
        otypen = 0;
    case 'modelout'
        otypen = 1;
    otherwise
        error('Unsupported output type');
end

if nargin < 5 && otypen ~= 0
    error('Coefs parameter has to be given to calculate output');
end

% check that x is column vector, fix if it is not
if size(x, 1) < size(x, 2)
    x = x.';
end


% precalculate abs(x)^k, x* and x^2 for faster implementation
abs_x = abs(x);
abs_x_powk = zeros(length(x), K);
temp = ones(size(abs_x));
x_conj = conj(x);
x2 = x.^2;

% go through all nonlinearity orders
% each iteration do multiplication (expecting it is faster than calculating
% the power)
abs_x_powk(:,1) = temp;
for k = 1:K-1
    temp = temp .* abs_x;
    abs_x_powk(:,k+1) = temp;
end


if otypen == 0
    coef_num = floor((K-1)/2)*(4*M+1)+(M+1);
    U = zeros(length(x),coef_num);
else
    U = zeros(size(x));
end
coef_ind = 1;

% matrix creation
for k = 0:((K-1)/2)
    for i = 0:M
        u = abs_x_powk(:,1+2*k) .* circshift(x,i);
        u(1:i) = 0;            % trimming values in history which are not available
        if otypen == 0
            U(:,coef_ind) = u;
        else
            U = U + coefs(coef_ind) .* u;
        end
        coef_ind = coef_ind + 1;
    end
end

for k = 1:((K-1)/2)
    for i = 1:M
        u = abs_x_powk(:,1+2*(k-1)) .* x2 .* circshift(x_conj,i);
        u(1:i) = 0;            % trimming values in history which are not available
        if otypen == 0
            U(:,coef_ind) = u;
        else
            U = U + coefs(coef_ind) .* u;
        end
        coef_ind = coef_ind + 1;
    end
end

% Simplified second order dynamics  
for k = 1:((K-1)/2)
    for i = 1:M
        u = abs_x_powk(:,1+2*(k-1)) .* x .* circshift(abs_x_powk(:,3),i);
        u(1:i) = 0;            % trimming values in history which are not available
        if otypen == 0
            U(:,coef_ind) = u;
        else
            U = U + coefs(coef_ind) .* u;
        end
        coef_ind = coef_ind + 1;
    end
end

for k = 1:((K-1)/2)
    for i = 1:M
        u = abs_x_powk(:,1+2*(k-1)) .* x_conj .* circshift(x2,i);
        u(1:i) = 0;            % trimming values in history which are not available
        if otypen == 0
            U(:,coef_ind) = u;
        else
            U = U + coefs(coef_ind) .* u;
        end
        coef_ind = coef_ind + 1;
    end
end
