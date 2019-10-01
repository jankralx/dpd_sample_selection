function [hplot] = PlotFrequencySpectrum(sig, hplot, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R = 50;

if nargin < 2
    hplot = 0;
end

global Fs;
if nargin < 3
    if exist('Fs', 'var')
        fs = Fs;
    else
        fs = 1;
    end
end

spect = fft(sig);
spect = 1/length(sig)*spect;
freq_axis = fftfreq(length(spect), fs, true);

if ~isa(hplot, 'matlab.graphics.chart.primitive.Line')
    hplot = plot(fftshift(freq_axis), db(fftshift(abs(spect)), R)+30);
else
    % +30 because of dBW -> dBm
    set(hplot, 'XData', fftshift(freq_axis), 'YData',db(fftshift(abs(spect)), R)+30);
end

