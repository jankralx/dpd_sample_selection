function [freqs] = fftfreq(N, fs, negative)
    if nargin < 3
        negative = 1;
    end
    
    freqs = ((0:N-1).')*fs/N;
    
    if negative
        freqs(1+floor((N+1)/2):end) = freqs(1+floor((N+1)/2):end)-fs;
    end    
end

