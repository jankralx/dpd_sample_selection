function [us] = Gradient_Sampling(A, b, n)
% [us] = GRADIENT_SAMPLING(A, X, N) returns sample selection vector
% based on the Gradient-based Sampling Algorithm of LS problems.
%
% The method was published here https://arxiv.org/pdf/1803.00841.pdf
%
%   Inputs:
%   =======
%   A  - matrix of equation coefficients
%
%   b - vector of right-hand-sides (A*x = b)
%
%   N  - number of samples to be selected
%
%   Returns:
%   ========
%
%   us - indexes of selected samples

% Authors: Jan Kral <kral.j@lit.cz>
% Date: 19.6.2019

% no selection if all samples need to be selected
if n >= size(A,1)
    us = 1:size(A,1);
    return
end

% calculate x0 (good guess) - it will be the full LS solution as we only
% need the reference results

% x0 = (A'*A)\(A'*b);
% x0 = lscov(A, b);
% x0 = zeros(size(A,2),1);

% we can optimize as x0 is zero vector

% calculate vector of |gradients g| (Eq. 4)
g = zeros(size(A,1),1);
for i = 1:size(A,1)
%     g(i) = norm(A(i,:)*(b(i)-A(i,:)*x0));
    g(i) = norm(A(i,:)*b(i));
end
% g = sum(abs(A).^2,2).*b;

% calculate probabilities p (Eq. 5)
p = g / sum(g);

% now select n samples randomly
us = zeros(size(A,1),1);
try
    sel = randsrc(n+10, 1, [1:length(p); p.']);      % randomly select indexes of elements with corresponding probabilities
    sel = unique(sel);
    if numel(sel) > n
        % randomly select indexes which will be excluded
        excl_ind = randperm(numel(sel), numel(sel)-n);
        sel(excl_ind) = [];
    end    

    us(sel) = 1;
    num_samples = numel(sel);
    i = 1;
    watchdog = 2;
    while num_samples < n && watchdog
        if ~us(i)           % current sample can be selected only if it is not already selected
            sel = rand < n*p(i);
            us(i) = sel;
            num_samples = num_samples + sel;
        end

        % move i to next row
        i = i + 1;
        if i > size(A,1)
            i = 1;
            watchdog = watchdog - 1;
        end
    end
catch
end
us = us > 0;
end
