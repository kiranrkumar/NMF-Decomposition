% STFT function to produce a result that will work with the provided
% iSpectrogram function
%
% Usage: [S, F, T] = myStft(x, nfft, overlap, winlen, fs);
%
% Parameters:
%   x:          input, mono, time-domain signal
%   nfft:       number of FFT bins to calculate
%   overlap:    stft overlap size in samples
%   winlen:     size of signal time frames in samples
%   fs:         sampling rate
%
% Returns:
%   S:          complex data spectrogram
%   F:          vector of calculated frequencies
%   T:          vector time increments in seconds

function [S, F, T] = myStft(x, nfft, overlap, winlen, fs)
    %make mono and normalize
    x = sum(x, 2) / size(x, 2);
    x = x / max(abs(x));
    
    %buffer the input audio and window
    xBuff = buffer(x, winlen, overlap, 'nodelay');
    winmat = repmat(hamming(winlen), [1, size(xBuff,2)]);
    xBuff = xBuff .* winmat;
    
    %take the FFT
    S = fft(xBuff, nfft);
    %T = 0 : (nfft-overlap)/fs : ((size(S, 2) - 1)*(nfft-overlap))/fs;
    
    %create [F]requency vector
    freqRes = fs/nfft;
    F = freqRes:freqRes:fs;
    
    %create [T]ime vector
    hopsizeSec = (winlen - overlap)/fs;
    totalTime = hopsizeSec * (size(S,2) - 1) + winlen/fs;
    T = 0:hopsizeSec:totalTime;
    
end