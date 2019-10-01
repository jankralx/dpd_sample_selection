function [ rx_bin, rx_symb_iq, rx_filt ] = QAM_Demodulator( rx_sig, params )
% QAM Demodulator
%   Universal QAM demodulator with optional Raised Cosine Filter
%   rx_sig - baseband received signal to be demodulated
%
%   Params:
%   =======
%
%   params.M - specifies order of QAM modulation (QAM16, QAM64 ...)
%              (default = 4)
%   params.filt = 'root' for root-raised cosine filter (default)
%               = 'normal' for raised cosine filter
%               = 'none' for no filter
%
%   params.rcf.rolloff - raised cosine roll-off factor (default = 0.15)
%   params.rcf.span    - raised cosine filter is truncated to span symbols
%                        (default = 80)
%
%   params.oversample_fact - oversampling factor (default = 8)
%
%   Returns:
%   ========
%
%   rx_bin     - demodulated binary data in vector
%   rx_symb_iq - sampled symbols in IQ format
%   rx_filt    - signal filtered by receiver filter (if any)

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


if ~strcmp(params.filt,'none')
    % RRC filter
    rc_filter = rcosdesign(params.rcf.rolloff, ...
        params.rcf.span, params.oversample_fact, params.filt);

    rx_filt = filter(rc_filter, 1, rx_sig);
    delay = length(rc_filter)-1;
else
    rx_filt = rx_sig;
    delay = 0;
end

% sample rx signal to obtain rx symbols
rx_symb_iq = rx_filt(1+delay:params.oversample_fact:end);

% denormalize
rx_symb_iq = rx_symb_iq * (sqrt(params.M)-1); 

rx_symb_dec = qamdemod(rx_symb_iq, params.M);
rx_symb_mat = de2bi(rx_symb_dec);  % convert to integers

% reshape demodulated data to obtain vector of bits
rx_bin = reshape(rx_symb_mat, numel(rx_symb_mat), 1);

end

