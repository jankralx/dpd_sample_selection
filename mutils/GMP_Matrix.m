function U = GMP_Matrix(x, Ka, La, Kb, Lb, Mb, Kc, Lc, Mc, type, otype, coefs)

% GMP_Matrix function returns matrix U created from the given signal x.
%   Matrix is created for Generalized Memory Polynomial (GMP) model
%   of nonlinearity order Ka, Kb, Kc, memory length La, Lb, Lc, and
%   number of cross-terms Mb, Mc.
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the DDR2 matrix U
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
%   U  - GMP matrix created from signal x (otype = 'matrix') or output
%        signal of GMP model with coefficients coefs (otype = 'modelout').


% Authors: Jan Kral <kral.j@lit.cz>
% Date: 10.2.2019

% the GMP is defined in the article references/morgan_GMP
% https://ieeexplore.ieee.org/abstract/document/1703853
% Equation 23 and 24

if nargin < 10
    if length(Ka) > 1 || length(La) > 1 || length(Kb) > 1 || ...
            length(Lb) > 1 || length(Mb) > 1 || length(Kc) > 1 || ...
            length(Lc) > 1 || length(Mc) > 1
        type = 'all';
    else
        type = 'num';
    end
end

if nargin < 11
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

if nargin < 12 && otypen ~= 0
    error('Coefs parameter has to be given to calculate output');
end

% check that x is column vector, fix if it is not
if size(x, 1) < size(x, 2)
    x = x.';
end

% convert GMP parameters to sets
if strcmp(type, 'all')
    ka = Ka;
    la = La;
    kb = Kb;
    lb = Lb;
    mb = Mb;
    kc = Kc;
    lc = Lc;
    mc = Mc;
else
    ka = 0:(Ka-1);
    la = 0:(La-1);
    kb = 1:Kb;
    lb = 0:(Lb-1);
    mb = 1:Mb;
    kc = 1:Kc;
    lc = 0:(Lc-1);
    mc = 1:Mc;    
end

% precalculate abs(x)^k for faster implementation
k_pre = union([ka, kb], kc);        % union of all non-linearity orders - all needed powers

abs_x = abs(x);
abs_x_powk = zeros(length(x), length(k_pre));
temp = ones(size(abs_x));

% go through all possible nonlinearity orders
% each iteration do multiplication (expecting it is faster than calculating
% the power) and save only those which are in k_pre - all needed powers
ind = 1;
for k = 0:max(k_pre)
    if ismember(k, k_pre)
        abs_x_powk(:,ind) = temp;
        ind = ind + 1;
    end
    temp = temp .* abs_x;
end

% prepare indexes to abs_x_powk for ka, kb, and kc
ka_pre = zeros(size(ka));
kb_pre = zeros(size(kb));
kc_pre = zeros(size(kc));
for i = 1:length(ka)
    ka_pre(i) = find(k_pre == ka(i));
end
for i = 1:length(kb)
    kb_pre(i) = find(k_pre == kb(i));
end
for i = 1:length(kc)
    kc_pre(i) = find(k_pre == kc(i));
end

% calculate number of coefficients
coef_num = length(ka)*length(la) + length(kb)*length(lb)*length(mb) + ...
    length(kc)*length(lc)*length(mc);

if otypen == 0
    U = zeros(length(x),coef_num);
else
    U = zeros(size(x));
end
coef_ind = 1;

% matrix creation
for k_ind = 1:length(ka)
    for lv = la
        u = circshift(x,lv).*circshift(abs_x_powk(:,ka_pre(k_ind)), lv);
        u(1:lv) = 0;            % trimming values in history which are not available
        if otypen == 0
            U(:,coef_ind) = u;
        else
            U = U + coefs(coef_ind) .* u;
        end
        coef_ind = coef_ind + 1;
    end
end

for k_ind = 1:length(kb)
    for lv = lb
        for m = mb
            u = circshift(x,lv).*circshift(abs_x_powk(:,kb_pre(k_ind)), lv+m);
            u(1:(lv+m)) = 0;    % trimming values in history which are not available
            if otypen == 0
                U(:,coef_ind) = u;
            else
                U = U + coefs(coef_ind) .* u;
            end
            coef_ind = coef_ind + 1;
        end
    end
end

for k_ind = 1:length(kc)
    for lv = lc
        for m = mc
            u = circshift(x,lv).*circshift(abs_x_powk(:,kb_pre(k_ind)), lv-m);
            u(1:lv) = 0;                % trimming values in history which are not available
            u(end+1+(lv-m):end) = 0;    % trimming values in future which are not available
            if otypen == 0
                U(:,coef_ind) = u;
            else
                U = U + coefs(coef_ind) .* u;
            end
            coef_ind = coef_ind + 1;
        end
    end
end
