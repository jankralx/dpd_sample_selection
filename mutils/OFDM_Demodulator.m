function [rx_bin, imod_symbs] = OFDM_Demodulator(rx_sig, par, ref_symbs)
% Universal demodulator for OFDM (and f-OFDM)
%
%   bin_in - binary data input
%
%   ref_symbs - reference symbols used by the internal demodulator, used to
%               adjust gain and phase rotation of demodulated OFDM symbols
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
%   rx_bin     - demodulated binary data
%   imod_symbs - demodulated symbols in IQ format (each row is one frame,
%                each column is one subcarrier)

if nargin < 2
    par = struct();
elseif nargin < 3
    if ~isa(par, 'struct')
        ref_symbs = par;
        par = struct();
    end
else
    if ~isa(par, 'struct')
        par = struct();
    end    
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


% calculate number of data subcarriers
num_data_carriers = par.num_rb * par.rb_size;

% make a column vector from rx_sig
rx_sig = rx_sig(:);

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

    %%% filter the input signal %%%
    rx_sig = filter(filt_num, 1, rx_sig);
    rx_sig = rx_sig(par.filt.L:end);           % trim filter delay
end

%% OFDM demodulator

% reshape to get frame signals
num_frames = length(rx_sig)/(par.fft_len+par.cp_len);

if mod(length(rx_sig),par.fft_len+par.cp_len)
    error('Length of input signal has to be multiple of frame length (par.fft_len+par.cp_len)');
end

ofdm_sig = reshape(rx_sig, par.fft_len+par.cp_len, num_frames).';

% remove cyclic prefix
fft_in = ofdm_sig(:,par.cp_len+1:end);

% fft to get OFDM symbols
ofdm_symbs = fftshift(fft(fft_in,[],2),2);

% select data subcarriers to get imod symbols
fft_offset = (par.fft_len - num_data_carriers) / 2; % for band center
imod_symbs = ofdm_symbs(:,fft_offset+(1:num_data_carriers));


%% internal demodulator

% set parameters for internal QAM modulator
qam_par.M = par.imod.M;         % set QAM - M state
qam_par.filt = 'none';          % no filter
qam_par.oversample_fact = 1;    % no oversampling

% reshape internal symbols to get sequence of symbols
imod_sym = reshape(imod_symbs, numel(imod_symbs), 1);

% currently there is no training sequence, so we use given reference to
% adjust gain and phase rotation
g = lscov(imod_sym, reshape(ref_symbs, numel(ref_symbs),1));

% correct gain for both calculated and received symbols
imod_symbs = g * imod_symbs;
imod_sym = g * imod_sym;

% demodulate bits from symbols
rx_bin = QAM_Demodulator(imod_sym, qam_par);



