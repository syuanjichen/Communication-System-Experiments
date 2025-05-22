clc; clear; close all;

%% OFDM MOD

% Define the FFT size
nfft = 64;
cp = 16;

% Given short training symbol in frequency domain
short_training_freq = sqrt(13/6)*[0, 0, 0, 0, -1-1j, 0, 0, 0, -1-1j, ...
                      0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1+1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, -1-1j, 0, 0, 0, ...
                      -1-1j, 0, 0, 0, 1+1j, 0, 0, 0];

%long_training_freq = [0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];
long_training_freq = [1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];

% Convert to time domain using IFFT
short_training_time = ifft(short_training_freq, nfft);
long_training_time = ifft(long_training_freq, nfft);

% Repeat the short training symbol (10 times, as in typical 802.11)
short_training_cp = short_training_time(nfft-cp+1:nfft);
long_training_cp = long_training_time(nfft-2*cp+1:nfft);

sts = [short_training_cp, short_training_time, short_training_cp, short_training_time].';
lts = [long_training_cp, long_training_time,long_training_time].';

FFT_size = 64;              % FFT size
cyclic_prefix_len = 16;     % Cyclic prefix length
fc = 900e6;                 % Carrier frequency in Hz
sampling_rate = 10e6;       % Sampling rate before upsampling (20 MHz)
mod_order = 4;              % Modulation order (16-QAM)
interpFactor = 1;
ppm = 5/1e6;


num_OFDM_symbols = 100; % Adjust this as desired
num_bits = FFT_size * log2(mod_order) * num_OFDM_symbols;

% 1. Generate random bits
rng(21);
bits = randi([0 1], num_bits, 1);

% 2. Modulate the bits with 16-QAM
dataSymbols = qammod(bits, mod_order, 'InputType', 'bit', 'UnitAveragePower', true);

% 3. Reshape and Map to OFDM symbols
dataSymbols = reshape(dataSymbols, FFT_size, num_OFDM_symbols);

% 4. OFDM Modulation
ofdmSymbols = ofdmmod(dataSymbols, FFT_size, cyclic_prefix_len);

% 3. Concatenate STS, LTS, and Data Symbols to Form Frame
frame = [sts; lts; ofdmSymbols];

% 5. Plot the Frame in Time Domain
time_axis = (0:length(frame)-1) * (1 / sampling_rate) * 1e6; % Convert to microseconds (µs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_vector = (0:length(frame)-1)/ sampling_rate ;

frame = frame .* exp(1j * 2 * pi * fc *(ppm)* time_vector).';

%% USRP define
txChannel = 'Channel 1';       % USRP transmit channel
rxChannel = 'Channel 2';       % USRP receive channel
Gain = 10;
radio_Tx = comm.SDRuTransmitter(...
                        'Platform',             "N200/N210/USRP2", ...
                        'IPAddress',            '192.168.10.2', ...
                        'CenterFrequency',      fc, ...
                        'Gain',                 Gain);
radio_Rx = comm.SDRuReceiver(...
                    'Platform',             "N200/N210/USRP2", ...
                    'IPAddress',            '192.168.10.2', ...
                    'CenterFrequency',      fc, ...
                    'Gain',                 Gain, ...
                    'OutputDataType',       'double' , ...
                    'SamplesPerFrame',  1.8*length(frame)  ...
                    );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USRP TX, RX
f = 1;
t = linspace(0, 3, 1000);
% data = 1*sin(2*pi*f*t).';
for ii = 1:20
    fprintf("Tx i = %03d", ii);
    tunderrun = radio_Tx(frame);
    [rcvdSignal, ~, toverflow] = step(radio_Rx);
end
%%
% figure(1);
% plot(real(rcvdSignal));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OFDM DeMod
fprintf("\nStart DeOfdm\n");
mfir = conj(flipud(lts));
timestep = (0:length(mfir)-1) * (1 / sampling_rate) * 1e6;

% filSig = conv(rcvdSignal, mfir, 'full');
filSig = filter(mfir, 1, rcvdSignal);

[m,idx] = max(abs(filSig));

fprintf("Start identify DeQAM\n");

%% CFO

slts = rcvdSignal(idx-nfft+1:idx);
flts = rcvdSignal(idx-2*nfft+1:idx-nfft);


p = sum(conj(flts) .* slts);

phi = atan2(imag(p),real(p));

f_ratio = phi / (2*pi*nfft); % delta_f / fs
delta_f = f_ratio * sampling_rate; 

our_ppm = delta_f / fc * 1e6;
fprintf("\n ppm = %f, \n our ppm = %f\n", ppm*1e6, our_ppm);
% raw no cfo correction
rx_sym_raw = ofdmdemod(rcvdSignal(idx+1:length(ofdmSymbols)+idx), FFT_size, cyclic_prefix_len);

newtime = [zeros(1, idx-320), time_vector, zeros(1, length(rcvdSignal)-length(time_vector)-idx+320)];

rcvdSignal = rcvdSignal.* exp(-1j * 2 * pi * delta_f .* newtime).';

%%
rxSymbols = ofdmdemod(rcvdSignal(idx+1:length(ofdmSymbols)+idx), FFT_size, cyclic_prefix_len);

rxSymbols_noeq = rxSymbols;
%% equalization
lts_rcvd = rcvdSignal(idx-nfft+1:idx).';
% h = long_training_freq .* fft(lts_rcvd, nfft) ;
h = fftshift(fft(lts_rcvd, nfft)) ./ fftshift(long_training_freq);

h = h.';

rxSymbols = rxSymbols ./ h;


rxBits = qamdemod(rxSymbols(:), mod_order, 'OutputType', 'bit', 'UnitAveragePower', true);
% Calculate Bit Error Rate (BER)
[numErrors, ber] = biterr(bits, rxBits);
fprintf('Bit Error Rate (BER): %f\n', ber);


%% QAM plot
% figure(3);
scatterplot(rxSymbols(:));
title('Received Constellation with CFO correction');
% figure(4);
scatterplot(rx_sym_raw(:));
title('Rx symbol without CFO correction');


scatterplot(rxSymbols_noeq(:));
title('Rx symbol with CFO correction, wihout Equalization');









