function [tx_sig, imod_symbs] = OFDM_Modulator(bin_in, par)
% Universal modulator for OFDM (and f-OFDM)
%
%   bin_in - binary data input
%
%   Params:
%   =======
%   par.fft_len - number of FFT points (default 1024)
%
%   par.num_rb - number of resource blocks (default 50)
%
%   par.rb_size - number of subcarriers per resource block (default 12)
%
%   par.cp_len - cyclic prefix length in samples (default 72)
%
%   par.filt.ena - if set true then f-OFDM modulator is used
%                  (default false)
%
%   par.filt.L - filter length (= filter order + 1), has to be odd
%                (default 513)
%
%   par.filt.tone_offset - tone offset or excess bandwidth (in subcarriers)
%                          (default 2.5)
%
%   par.imod.name - specifies internal modulation ('QAM')
%                   (default 'QAM')
%
%   par.imod.name == 'QAM':
%   par.imod.M - specifies order of internal QAM modulation (QAM16, QAM64 ...)
%                (default = 64)
%   
%
%   Returns:
%   ========
%
%   tx_sig     - modulated signal
%   imod_symbs - modulated symbols in IQ format (each row is one frame,
%                each column is one subcarrier)

if nargin < 2
    par = struct();
end

% register input parameters
par = RegisterParam(par, 'fft_len', 1024);
par = RegisterParam(par, 'num_rb', 50);
par = RegisterParam(par, 'rb_size', 12);
par = RegisterParam(par, 'cp_len', 72);
par = RegisterParam(par, 'filt');
par.filt = RegisterParam(par.filt, 'ena', false);
par.filt = RegisterParam(par.filt, 'L', 513);
par.filt = RegisterParam(par.filt, 'tone_offset', 2.5);

par = RegisterParam(par, 'imod');
par.imod = RegisterParam(par.imod, 'name', 'QAM');
par.imod = RegisterParam(par.imod, 'M', 64);

% for F-OFDM check filter length is fft_len/2 + 1 - otherwise it will lead
% to suboptimal EVM
if par.filt.ena && par.filt.L ~= par.fft_len / 2 + 1
    warning('par.filt.L != fft_len/2 + 1 which leads to suboptimal EVM');
end

% calculate number of data subcarriers
num_data_carriers = par.num_rb * par.rb_size;


%% modulate input binary data to QAM symbols

% set parameters for internal QAM modulator
qam_par.M = par.imod.M;         % set QAM - M state
qam_par.filt = 'none';          % no filter
qam_par.oversample_fact = 1;    % no oversampling

bin_in = bin_in(:);                                     % convert input bits into vector
bits_per_frame = num_data_carriers * log2(qam_par.M);   % calculate number of bits per frame
num_frames = ceil(length(bin_in)/bits_per_frame);       % number of frames needed to transmitt all input bits
num_align_bits = num_frames * bits_per_frame - length(bin_in);

if num_align_bits ~= 0
    % needs to add zero bits for alignment
    bin_in = [bin_in; zeros(num_align_bits,1)];
    warning('Adding alignemt bits into sent binary data -> Set length of bin_in to be multiple of par.num_rb * par.rb_size * log2(par.imod.M)');
end

% get QAM symbols
[~, imod_symbs] = QAM_Modulator(bin_in, qam_par);
% reshape symbols such each row is one frame and each column is one subcarrier
imod_symbs = reshape(imod_symbs, num_frames, num_data_carriers);    

%% OFDM modulator

fft_offset = (par.fft_len - num_data_carriers) / 2; % for band center
ofdm_symbs = [zeros(num_frames, fft_offset) imod_symbs ...
              zeros(num_frames, par.fft_len-fft_offset-num_data_carriers)];

ifft_out = ifft(ifftshift(ofdm_symbs,2),[],2);

% add cyclic prefix
ofdm_sig = [ifft_out(:,end-par.cp_len+1:end) ifft_out];

% reshape to get linear signal
tx_sig = reshape(ofdm_sig.', numel(ofdm_sig), 1);

%% f-OFDM part

if par.filt.ena
    %%% design prototype filter %%%
    % as implemented here https://www.mathworks.com/help/comm/examples/f-ofdm-vs-ofdm-modulation.html#OFDMvsFOFDMExample-11
    % ref: R1-165425. "F-OFDM scheme and filter design." 3GPP TSG RAN WG1 meeting 85. Huawei; HiSilicon. May 2016.
    %
    half_filt = floor(par.filt.L/2);
    n = -half_filt:half_filt;

    % sinc prototype filter
    pb = sinc((num_data_carriers+2*par.filt.tone_offset) .* n ./ par.fft_len);

    % sinc truncation window
    w = (0.5*(1+cos(2*pi.*n/(par.filt.L-1)))).^0.6;

    % normalize filter coefficients
    filt_num = (pb.*w) / sum(pb.*w);

    %%% filter the OFDM signal %%%
    % filter with zero-padding to flush tail
    tx_sig = filter(filt_num, 1, [tx_sig; zeros(par.filt.L-1,1)]);
end




