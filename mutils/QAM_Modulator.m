function [tx_sig, tx_symb_iq] = QAM_Modulator(bin_in, params)
% QAM modulator
%   Universal QAM modulator with optional Raised Cosine Filter
%   bin_in - binary data input
%
%   Params:
%   =======
%
%   params.M - specifies order of QAM modulation (QAM16, QAM64 ...)
%              (default = 4)
%   params.filt = 'sqrt' for root-raised cosine filter (default)
%               = 'normal' for raised cosine filter
%               = 'none' for no filter
%
%   params.rcf.rolloff - raised cosine roll-off factor (default = 0.15)
%   params.rcf.span    - raised cosine filter is truncated to span symbols
%                        (default = 80)
%
%   params.oversample_fact - oversampling factor (default = 8)
%                            oversampling work only if RC filter enabled
%
%   Returns:
%   ========
%
%   tx_sig     - modulated signal
%   tx_symb_iq - modulated symbols in IQ format

if nargin < 2
    params = struct();
end

% register input parameters
params = RegisterParam(params, 'M', 4);
params = RegisterParam(params, 'filt', 'sqrt');
params = RegisterParam(params, 'rcf');
params.rcf = RegisterParam(params.rcf, 'rolloff', 0.15);
params.rcf = RegisterParam(params.rcf, 'span', 80);
params = RegisterParam(params, 'oversample_fact', 8);

k = log2(params.M);                % number of bits per symbol

% reshape data to form data symbols
data_in_mat = reshape(bin_in, length(bin_in)/k, k);
tx_symb_dec = bi2de(data_in_mat);  % convert to integers

tx_symb_iq = qammod(tx_symb_dec, params.M);

% normalize such that maximum signal value in real or imaginary axis is 1
tx_symb_iq = tx_symb_iq ./ (sqrt(params.M)-1); 

if ~strcmp(params.filt,'none')
    % RRC filter
    rc_filter = rcosdesign(params.rcf.rolloff, ...
        params.rcf.span, params.oversample_fact, params.filt);

    % upsample signal
    tx_sig = upfirdn(tx_symb_iq, rc_filter, params.oversample_fact, 1);
else
    tx_sig = tx_symb_iq;
end

end
