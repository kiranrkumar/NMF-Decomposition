% Note: This function works on an STFT that DOES contain the redundant half
% of the complex spectrum.
% iSpectrogram Performs the Inverse Short Time FFT
%                                                                                                                         
% Usage: [data] = iSpectrogram(B, nfft, winlen, overlap);
% 
% Parameters:
%   B:          [nfft x timeSegments] STFT data
%   NFFT:       Number of points for which the FFT was calculated per
%               slice
%   fs:         Sampling frequency of the original data
%   winlen:     window length in samples
%   overlap:    number of samples for whic consecutive FFT slices overlap
%   [winvect]:  nfft x 1 windowing function vector that was used to window 
%               the time segments of the signal before FFT'ings
%   [cmpxData]: the complex STFT data from the original time domain signal. 
%               Required if B is a magnitude spectrum without complex data.
%
% Returns:
%   data:     column vector containing the time domain representation of
%               the input STFT B


function [data] = iSpectrogram(B, nfft, winlen, overlap, winvect, cmpxData)

    if (nargin == 5 && isreal(B))
        error('You must provide phase data if processing a magnitude spectrum');
    elseif (nargin == 6 && isreal(cmpxData))
        error('cmpxData must be complex data');
    %incorporate phase data into the magnitude spectrum. B must be real
    %(magnitude spectrum), and cmpxData must be complex numbers (FFT results)
    elseif (nargin == 6 && ~isreal(cmpxData) && isreal(B))

        %Calculate the desired complex numbers by scaling each original complex
        %number by the newly calculated magnitude value at that position
        B = (B./abs(cmpxData)) .* cmpxData;

        %A more complicated but valid way of calculating phase
        %     reals = real(cmpxData);
        %     imags = imag(cmpxData);
        %     angles = atan(imags ./ reals);
        %     
        %     % adjust angles as necessary since atan only gives -90 <= x <= 90
        %     for m = 1:size(angles,1)
        %         for n = 1:size(angles,2)
        %             if (reals(m,n) < 0)
        %                 angles(m,n) = angles(m,n) + pi;
        %             end
        %         end
        %     end
        %     
        %     B = B .* cos(angles) + B .* sin(angles) * 1i;

    end

    hopsize = winlen - overlap;
    [a,b]=size(B);

    % create a vector to hold the single column time-domain signal
    ispecgram = zeros((((hopsize*(b-1))+a)),1);

    B = ifft(B,nfft);

    %unwindow the signal if a window function was specified
    if nargin == 5
        winmat = repmat(winvect, [1, b]);
        B = B ./ winmat;
    end

    ispecStart = 1;
    ispecEnd = ispecStart - 1 + winlen;
    for curTimeCol = 1:1:b
        ispecgram(ispecStart:ispecEnd,1) = B(1:winlen, curTimeCol) + ispecgram(ispecStart:ispecEnd,1);
        ispecStart = ispecStart + hopsize;
        ispecEnd = ispecStart - 1 + winlen;
    end
    NX =(b*(winlen-overlap)) + overlap;
    data =  real(ispecgram);
    data = data(1:NX);

end

