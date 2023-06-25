N = 16;
halfN = floor((N+1)/2);
% Divide a single period among N segments, starting at 0
t = 2*pi*(0:(N-1))/N;

% Create sine samples at frequencies 1 to halfN, and cosine samples at
% frequencies 0 to halfN
sines = sin((1:(halfN-1))'*t);
cosines = cos((0:halfN)'*t);

% Combine the sines and cosines into a basis matrix
basis = [sines; cosines];

% Create sines of 3 different frequencies
sin1 = sin(t+pi/4);
sin2 = sin(2*t-pi/2);
sin3 = sin(3*t+pi);

% Combine the sines into a single vector
f = sin1+sin2+sin3;
F = f*basis';

% Calculate the amplitude and phase based on the sine and cosine weights
% from F.

Amp = sqrt(F(1:(halfN-1)).^2+F(halfN+1:N-1).^2);
Phas = atan2(F(halfN+1:N-1),F(1:(halfN-1)));

%

% Since              
%             F = f*basis'
% F*inv(basis') = f*basis'*inv(basis')
% F*inv(basis') = f
%
% The column vectors of basis are orthogonal, which guarantees that it is
% invertible.

fnew = F/basis';

%

fFft = fft(f);
ampFft = fftshift(abs(fFft));
phaseFft = fftshift(angle(fFft));

% The signs of the phases have some differences, but they have the correct
% magnitudes. Similarly, the amplitudes are correct, although everything is
% shifted so that the 0 is at the center.

%

iSftTest = fft(1:N);
isft(iSftTest);
isftMat(iSftTest);

%

t = 2*pi*(0:15)/16;
a = rand(1,3);
p = 2*pi*rand(1,3)-pi;
x = a(1)*cos(t+p(1))+a(2)*cos(3*t+p(2))+a(3)*cos(5*t+p(3));

[ampX,phaseX] = interpFft(fft(x));

%

dispSpectrum(x,16)

%

randX = rand(1,64);
symX = [0 randX(1:31) 0 flip(randX(1:31))];
t64 = (0:63);

figure()
subplot(2,2,1)
plot(t64,fft(symX))
title("Symmetric")

% There is a single peak at 0 and the rest of the values are low, with some
% fluctuation.

asymX = [0 randX(1:31) 0 -flip(randX(1:31))];
subplot(2,2,2)
plot(t64,fft(asymX))
title("Antisymmetric")

% The antisymmetric signal is almost all 0s, since there are no cosines in
% our function.

subplot(2,2,3)
plot(t64,fft(randX))
title("Random")

% The random signal behaves similarly to the symmetric one, with a single
% high point at 0 and low, fluctuating values for all other frequencies.

randXi = rand(1,64,"like",i);
subplot(2,2,4)
plot(t64,randXi)
title("Imaginary")


fMag = abs(N/2-(1:N));
fPhase = 2*pi*rand(1,N/2)-pi;
fPhase = [fPhase flip(fPhase)];

Freal = fMag.*cos(fPhase);
Fimag = fMag.*sin(fPhase);

F2 = Freal+i*Fimag;

f2 = ifft(F2);

fftMag = abs(fft(f2));

fShift = circshift(f2,N/2);
fShiftFft = fft(fShift);
fShiftMag = abs(fShiftFft);
fShiftPhase = angle(fShiftFft);

fShift2 = circshift(f2,N/4);
fShiftFft2 = fft(fShift2);
fShiftMag2 = abs(fShiftFft2);
fShiftPhase2 = angle(fShiftFft2);

fShiftPhase3 = angle(fft(circshift(f2,N/8)));

% The magnitudes are the same for all circular shifts. The phases are the
% same for the first frequency, and overlap for some other values, but
% other than that are different and don't have any distinguishable pattern.


fRev = flip(f2);
fRevFft = fft(fRev);
fRevMag = abs(fRevFft);
fRevPhase = angle(fRevFft);

% Reversing f2 similarly doesn't show any discernible pattern, except for
% overlap for the first value.


% g[n] = f2[n]*(-1)^n
g = f2.*((-1).^(1:N));
G = fft(g);
gMag = abs(G);
gPhase = angle(G);

t16 = (0:(N-1));
figure(); hold on
plot(t16,fPhase)
plot(t16,gPhase)
title("f2 vs. g"); hold off

% Looking at the vectors fMag and gMag, gMag appears to be equal to fMag,
% shifted circularly by N/2. Graphically, the phases of f2 and g appear to
% be reflections across the line y=4. 
%
% Again, the phase doesn't appear to have a pattern, although it looks
% roughly symmetrical for both.


% ======================================================================= %
%                               FUNCTIONS                                 %
% ======================================================================= %

function X = sft(x)
    N = length(x);
    X = 1:N; % just to create a vector of the correct size
    t = 2*pi*(0:(N-1))/N;
    for k = 0:(N-1)
      X(k+1) = x*exp(-i*k*t).';
    end
end


function [X,tMat] = sftMat(x)
    N = length(x);
    t = 2*pi*(0:(N-1))/N;
    k = (0:(N-1))';
    tMat = exp(-i*k*t);
    X = x*tMat;
end


function x = isft(X)
    N = length(X);
    x = 1:N;
    t = 2*pi*(0:(N-1))/N;
    for k = 0:(N-1)
        x(k+1) = X*exp(i*k*t).'/N;
    end
end

function x = isftMat(X)
    N = length(X);
    t = 2*pi*(0:(N-1))/N;
    k = (0:(N-1))';
    tMat = exp(i*k*t)/N;
    x = X*tMat;
end


function [amp,phase] = interpFft(X)
    amp = sqrt(imag(X).^2+real(X).^2);
    phase = atan2(imag(X),real(X));
end


function dispSpectrum(x,Fs)
    N = length(x);
    N2 = ceil(N/2);
    Nyquist = ceil(Fs/2);
    xFft = fft(x);

    amp = abs(xFft(1:N2));
    phase = angle(xFft(1:N2));
    freq = Nyquist*(0:(N2-1))/N2;

    figure()
    subplot(2,1,1); plot(freq,amp)
    axis("auto"); title("Amplitude"); xlabel("Frequency")
    subplot(2,1,2); plot(freq,phase)
    axis("auto"); title("Phase"); xlabel("Frequency")
end