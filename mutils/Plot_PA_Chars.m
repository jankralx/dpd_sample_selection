function [hpl_am, hpl_pm] = Plot_PA_Chars(pa_in, pa_out, hpl_am, hpl_pm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    hpl_am = 0;
    hpl_pm = 0;
end

% correct the phase offset
% pa_out = pa_out * exp(1j*angle(lscov(pa_out, pa_in)));

% calculate pm characteristics
ph_diff = angle(pa_out) - angle(pa_in);
ph_diff(ph_diff > pi) = ph_diff(ph_diff > pi) - 2*pi;
ph_diff(ph_diff < -pi) = ph_diff(ph_diff < -pi) + 2*pi;


% AM-PM
if ~isa(hpl_pm, 'matlab.graphics.chart.primitive.Line')
    subplot(2,1,2);
    hpl_pm = plot(abs(pa_in), ph_diff, '.');
else
    set(hpl_pm, 'XData', abs(pa_in), 'YData', ph_diff);
end

% AM-AM
if ~isa(hpl_am, 'matlab.graphics.chart.primitive.Line')
    subplot(2,1,1);
    hpl_am = plot(abs(pa_in), abs(pa_out), '.');
else
    set(hpl_am, 'XData', abs(pa_in), 'YData', abs(pa_out));
end
