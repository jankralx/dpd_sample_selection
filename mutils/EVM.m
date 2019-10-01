function [evm] = EVM(x, ref, dim, evm_type, norm_type)
% EVM calculates Error Vector Magnitude for given signal x using reference
% signal ref.
%  x   - symbols demodulated at receiver, vector or matrix, or
%        mutlidimensional array
%
%  ref - reference symbols, needs to be the same dimensions as x
%
%  dim - dimension across to calculate EVM, 0 to calculate EVM across all
%        dimensions (avarage of all samples)
%
% evm_type - ['rms' | 'peak']
%          - 'rms' calculates EVM using RMS of error vector magnitudes
%          - 'peak' calculates EVM using the maximum error vector magnitude
%          (default 'rms')
%
% norm_type - ['rms' | 'peak'] type of EVM normalization
%           - 'rms' normalize EVM using RMS of the reference symbols
%           - 'peak' normalize EVM using maximum vector of ref symbols
%           (default 'rms')

% Authors: Jan Kral <jankralx@gmail.com>
% Date: 15.3.2019

if nargin < 3
    dim = 0;
end
if nargin < 4
    evm_type = 'rms';
end
if nargin < 5
    norm_type = 'rms';
end

if dim == 0
    % convert both input signal and reference into vectors
    x = x(:);
    ref = ref(:);
    dim = 1;
end


if strcmp(evm_type, 'peak')
    err = sqrt(max(abs(ref-x).^2, [], dim));
else
    err = sqrt(mean(abs(ref-x).^2, dim));
end

if strcmp(norm_type, 'peak')
    nfact = sqrt(max(abs(ref).^2, [], dim));
else
    nfact = sqrt(mean(abs(ref).^2, dim));
end

% calculate EVM
evm = err ./ nfact;

