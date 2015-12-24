%%Decomposition and Reconstruction of an audio signal via non-negative
%matrix factorization.
%
%[W,H] = nmf_audioDecomp() returns the estimated basis and activation
%matrices of a default audio signal using MATLAB's built-in nnmf function
%(with k = 30).
%
%[W,H] = nmf_audioDecomp(audiofile, kOverride) allows the user to specify
%their own audio file to analyze. kOverride provides a method to override
%the default k value of 30 basis functions.
%
%[W, H, MAG, resultTimeDom, residual] = nmf_audioDecomp() returns the
%magnitude spectrogram MAG of the computed matrix W * H. resultTimeDom
%contains the time domain version of the computed W*H signal, and residual
%contains the difference of the original (signal - resultTimeDom).
%
%% Example 1:
%   To retrieve data on the default audio signal
%   'PapaWasARollingStoneMIDI.aif'
%
%   [W, H, MAG, resultTimeDom, residual] = nmf_audioDecomp;
%
%% Example 2:
%   To analyze an audio file called 'mySong.wav' and search for a maximum
%   of 50 components.
%
%   filename = 'mySong.wav'                   % filepath to the audio
%   k = 50;                                   % 50 components
%   [W, H] = nmf_audioDecomp(filename, k);    % Find basis and activations
%
%See also NNMF
function [W, H, MAG, resultTimeDom, residual] = nmf_audioDecomp(audiofile, kOverride)
    close all; %close all existing windows
    
    if nargin == 0 %work with the default audio signal
        root = pwd;
        filename = 'PapaWasARollingStoneMIDI.aif';
        [x, fs] = audioread(fullfile(root, filename));
    else %else use a user specified signal
        [x, fs] = audioread(audiofile);
    end
    
    %max out at 30 seconds
    dur = 30;
    t = (1:dur*fs)';
    if size(x,1) < dur*fs
        t = t(1:size(x,1));
    end
    xSeg = x(t);
    t = t / fs;
    t = t - 1/fs;
    
    %force to mono by averaging all columns into 1 channel
    xSeg = mean(xSeg, 2);
    xSeg = xSeg / max(abs(xSeg));

    %% Get and plot spectrogram data

    winLength = 1024;
    overlapLength = 128;
    nfft = winLength;
    % get the spectrogram data
    [S, F, T] = myStft(xSeg, nfft, overlapLength, winLength, fs);
    %create a local hamming window since myStft uses a hamming window
    winvect = hamming(winLength);
    %magnitude spectrum
    MAG = abs(S);
    
    %% Do the NMF

    %K value
    K = 30;
    if nargin == 2
        K = kOverride;
    end
    [W, H] = nnmf(MAG, K);
    KVect = 1:K;

    %% Plot

    %NMF data
    f = figure('name', 'NMF Data');
    subplot(2, 2, 1);
    imagesc(T, F, mag2db(MAG));
    set(gca, 'YDir', 'normal');
    title('Original Audio');
    xlabel(['Time (0 - ', num2str(T(end)), ' seconds)']);
    ylabel(['Freq (0 - ', num2str(F(end)), ' Hz)']);

    subplot(2, 2, 2);
    imagesc(T, F, mag2db(W*H));
    title('Calculated W * H');
    set(gca, 'YDir', 'normal');
    xlabel(['Time (0 - ', num2str(T(end)), ' seconds)']);
    ylabel(['Freq (0 - ', num2str(F(end)), ' Hz)']);

    subplot(2, 2, 3);
    imagesc(KVect, F, W);
    title('W - Basis Functions');
    set(gca, 'YDir', 'normal');
    xlabel(['Signal Components: ', num2str(K),]);
    ylabel('Frequency Content');

    subplot(2, 2, 4);
    imagesc(T, KVect, H);
    title('H - Activations');
    set(gca, 'YDir', 'normal');
    xlabel('Activation Time (sec)');
    ylabel('Signal Components');

    pos = f.Position;

    %resize and window and change position
    pos(1) = pos(1)/5; %bring closer to left
    pos(3) = pos(3) * 2.2; %increase width
    pos(4) = pos(4) * 1.7; %increase height
    f.Position = pos;

    f2 = figure('name', 'Pre|Post Inverse Spectrogram');
    subplot(211);
    plot(t, xSeg);
    title('Original Signal');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    Sinv = iSpectrogram(S, nfft, winLength, overlapLength, winvect);
    subplot(212);
    Sinv = Sinv(1:length(t)); %set to the same length as the original
    plot(t, Sinv);
    title('Reconstructed Signal');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    f2.Position = pos;


    colormap parula;

    %Return the reconstructed result and the residual
    resultTimeDom = iSpectrogram((W*H), nfft, winLength, overlapLength, winvect, S);
    residual = iSpectrogram(MAG - W*H, nfft, winLength, overlapLength, winvect, S);

    %% Analysis
    resultTimeDom = resultTimeDom(1:length(t));

    %normalize audio
    resultTimeDom = resultTimeDom / max(abs(resultTimeDom));
    xSeg = xSeg / max(abs(xSeg));

    %Frequency Domain
    normMAG = norm(MAG);
    normWH = norm(W*H);
    magRes = MAG - W*H;
    normMR = norm(magRes);

    %Time Domain
    normX = norm(xSeg);
    normTimeCalc = norm(resultTimeDom);
    timeRes = xSeg - resultTimeDom;
    normTR = norm(timeRes);

    fprintf('\nSpectrogram Data\n=========================\n');
    fprintf('Magnitude norm: %f\nW*H norm: %f\nResidual Norm: %f\n\n', normMAG, normWH, normMR);
    fprintf('Time Domain Data\n=========================\n');
    fprintf('Original Signal: %f\nCalculated Signal: %f\nResidual Norm: %f\n', normX, normTimeCalc, normTR);
end