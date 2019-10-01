function [ spect, freq_axis ] = FrequencySpectrum( signal, Fs_in, R_load )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global Fs;

if nargin > 1
    fs = Fs_in;    
else
    fs = Fs;
end

if nargin < 3
    R_load = 50;
end

spect = 1/length(signal)*fft(signal);
freq_axis = fftfreq(length(spect), fs, true);

% if ~isa(obj.spect_trace, 'matlab.graphics.chart.primitive.Line')
    figure();
    obj.spect_trace = plot(fftshift(freq_axis), db(fftshift(abs(spect)), R_load)+30);
    title('Frequency Spectrum');
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dBm)');
% else
%     % +30 because of dBW -> dBm
%     set(obj.spect_trace, 'XData',fftshift(freq_axis), 'YData',db(fftshift(abs(Signal)), obj.R)+30);
% end


end


