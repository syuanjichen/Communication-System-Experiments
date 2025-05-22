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
frame = [zeros(100, 1); sts; lts; ofdmSymbols; zeros(100, 1)];

% 5. Plot the Frame in Time Domain
time_axis = (0:length(frame)-1) * (1 / sampling_rate) * 1e6; % Convert to microseconds (µs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
figure(1);
plot(real(rcvdSignal));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OFDM DeMod
fprintf("\nStart DeOfdm\n");
mfir = conj(flipud(lts));
timestep = (0:length(mfir)-1) * (1 / sampling_rate) * 1e6;

% filSig = conv(rcvdSignal, mfir, 'full');
filSig = filter(mfir, 1, rcvdSignal);

[m,idx] = max(abs(filSig));

fprintf("Start identify DeQAM\n");

%%
rxSymbols = ofdmdemod(rcvdSignal(idx+1:length(ofdmSymbols)+idx), FFT_size, cyclic_prefix_len);

rxSymbols_noeq = rxSymbols;
%% equalization
lts_rcvd = rcvdSignal(idx-nfft+1:idx).';
% h = long_training_freq .* fft(lts_rcvd, nfft) ;
h = fftshift(fft(lts_rcvd, nfft)) ./ fftshift(long_training_freq);

h = h.';
nvar = var(rcvdSignal(idx+1:idx+nfft) - mean(rcvdSignal(idx+1:idx+nfft)));
% nvar = repmat(nvar, nfft, num_OFDM_symbols);
heff = repmat(h, 1, num_OFDM_symbols);

% rxSymbols = ofdmdemod(rcvdSignal(idx+1:length(ofdmSymbols)+idx), FFT_size, cyclic_prefix_len);
rxSymbols = rxSymbols ./ h;
% rxSymbols = rxSymbols ./ (h + nvar ./ abs(h).^2);
% rxSymbols = ofdmEqualize(rxSymbols, heff, nvar);


rxBits = qamdemod(rxSymbols(:), mod_order, 'OutputType', 'bit', 'UnitAveragePower', true);
% Calculate Bit Error Rate (BER)
[numErrors, ber] = biterr(bits, rxBits);
fprintf('Bit Error Rate (BER): %f\n', ber);


%% QAM plot
% figure(3);
scatterplot(rxSymbols(:));
title('Received Constellation');
% figure(4);
scatterplot(rxSymbols_noeq(:));
title('No equalization Constellation');
%%
figure(5);
plot(real(rcvdSignal)/ max(abs(rcvdSignal)));
hold on;
plot(abs(filSig)/m);
% plot((idx+1 -320 : idx + 640), real(frame(101:end-100)));
hold off;
legend;
%% plot for lab report
figure(10);

% Time vector for x-axis in microseconds
time_axis = (0:length(rcvdSignal)-1) * (1 / sampling_rate) * 1e6; % Convert to microseconds (µs)

% Plot the real part of the received signal
plot(time_axis, real(rcvdSignal), 'b');
hold on;

% Add vertical lines with labels
xline((idx - 320) * (1 / sampling_rate) * 1e6, '--r', 'STS'); % Convert index to time in µs
xline((idx - 160) * (1 / sampling_rate) * 1e6, '--g', 'LTS'); % Convert index to time in µs
xline(idx * (1 / sampling_rate) * 1e6, '--k', 'Data start');  % Convert index to time in µs
xline((idx + 400) * (1 / sampling_rate) * 1e6, '--m', 'Data end');    % Adjust label if needed

% Plot settings
xlabel('Time (\muS)');
ylabel('Real Part');
title('Fig. a: Generated OFDM Symbol with Cyclic Prefix (Before Upsampling)');
grid on;
hold off;
%
% figure(11);









