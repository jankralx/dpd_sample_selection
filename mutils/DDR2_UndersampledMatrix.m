function U = DDR2_UndersampledMatrix(x, P, M, us)

% DDR2_matrix function returns undersampled matrix U created from the given signal x.
%   Matrix is created for Dynamic Deviation Reduction (DDR2) model of nonlinearity order P
%   and memory length M. Matrix is created for undersampled DDR2
%   calculation with given undersampling factor UsF
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the DDR2 matrix U
%
%   P  - DPD nonlinearity order
%
%   M  - DPD memory length
%
%   [us] - undersampling:
%        scalar as undersampling factor - number of samples of signal x which are
%         omitted in matrix (signal time shift is not influenced)
%        vector as undersampling matrix - non-zero elements define which samples of x signals are
%        taken into resulting DDR2 matrix (signal time shift is not influenced)
%
%   Returns:
%   ========
%
%   U  - DDR2 matrix created from signal x


% Authors: Jan Kral <kral.j@lit.cz>, Tomas Gotthans
% Date: 22.1.2018


% the DDR2 is defined in the article
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5773325 
% Equation 3

if nargin < 4
    us = 1;         % if undersampling is not given take all samples
end

x_shifted = zeros(length(x),M);
for q = 0:M
    x_shifted(:,q+1) = circshift(x, q);
    x_shifted(1:q,q+1) = 0;
end

if length(us) == 1
    us = 1:us:length(x);
end

% undersample the x signal and x_shifted matrix
x_shifted = x_shifted(us,:);
x = x(us);

coef_num = floor((P-1)/2)*(4*M+1)+(M+1);
U = zeros(length(x),coef_num);
coef_ind = 1;

% separated for better understanding
for k = 0:((P-1)/2)
    for i = 0:M
        U(:,coef_ind) = (abs(x).^(2*k)).*(x_shifted(:,i+1));
        coef_ind = coef_ind + 1;
    end
end

for k = 1:((P-1)/2)
    for i = 1:M
        U(:,coef_ind) = (abs(x).^(2*(k-1))).*(x.^2).*conj(x_shifted(:,i+1));
        coef_ind = coef_ind + 1;
    end
end

% Simplified second order dynamics  
for k = 1:((P-1)/2)
    for i = 1:M
        U(:,coef_ind) = (abs(x).^(2*(k-1))).*(x).*abs(x_shifted(:,i+1)).^2;
        coef_ind = coef_ind + 1;
    end
end

for k = 1:((P-1)/2)
    for i = 1:M
        U(:,coef_ind) = (abs(x).^(2*(k-1))).*conj(x).*x_shifted(:,i+1).^2;
        coef_ind = coef_ind + 1;
    end
end
