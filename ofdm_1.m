clear all; 
close all; 
clc;

%Initialization
numFFT = 1024;           % Number of FFT points
numRB = 50;              % Number of resource blocks
rbSize = 12;             % Number of subcarriers per resource block
cp = 72;                 % Cyclic prefix length in samples
bitsPerSubCarr = 4;      % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 18;              % SNR in dB
toneOffset = 2.5;        % Tone offset or excess bandwidth (in subcarriers)
L = 513;                 % Filter length (=filterOrder+1), odd

numDataCarr = numRB*rbSize;    % number of data subcarriers in subband
halfFilt = floor(L/2);
n = -halfFilt:halfFilt;
pb = sinc((numDataCarr+2*toneOffset).*n./numFFT);   %Sinc function prototype filter
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;  % Sinc truncation window
fnum = (pb.*w)/sum(pb.*w);      % Normalized lowpass filter coefficients
filtTx=fnum;

% Modulation
bitsIn = randi([0 1], bitsPerSubCarr*numDataCarr, 1);   % Generate data symbols
symbolsIn = qammod_1(bitsIn,2^bitsPerSubCarr);
offset = (numFFT-numDataCarr)/2;
symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                 zeros(numFFT-offset-numDataCarr,1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));       % Pack data into an OFDM symbol

txSigOFDM = [ifftOut(end-cp+1:end); ifftOut];   % Prepend cyclic prefix
txSigFOFDM = conv2(filtTx,txSigOFDM);   % Filter;
[M,N]=size(txSigFOFDM);
txSigFOFDM=reshape(txSigFOFDM,1,M*N);

% Through channel
rxSig = filter([1 -1 0.4],[-2.4 3.6],txSigFOFDM);
rxSig = awgn(rxSig, snrdB, 'measured');     % Add AWGN
filtRx=fliplr(fnum);                    % Generate receiver filter
filteredRx=conv2(filtRx,txSigFOFDM);

% Demodulation
rxSymbol = filteredRx(cp+1:end);        % Remove cyclic prefix
RxSymbols = fftshift(fft(rxSymbol));    % Perform FFT
dataRxSymbols = RxSymbols(offset+(1:2*numDataCarr));    % Select data subcarriers
demodRx=qamdemod_1(dataRxSymbols,2^bitsPerSubCarr);

scatterplot(bitsIn,2^bitsPerSubCarr)
scatterplot(txSigFOFDM,2^bitsPerSubCarr)
scatterplot(demodRx,2^bitsPerSubCarr)

% Filter impulse response
h = fvtool(fnum, 'Analysis', 'impulse', ...
    'NormalizedFrequency', 'off', 'Fs', 15.36e6); 

% Plot power spectral density (PSD)
[psd,f] = periodogram(txSigFOFDM, rectwin(length(txSigFOFDM)), ...
                      numFFT*2, 1, 'centered');
plot(f,10*log10(psd));

%Bit Error Rate
%[noe,ber] = biterr(bitsIn,demodRx);
