clear all;
clc;

fft_size = 64;
cyclic_prefix_size = 16;
carrier_freq = 2.4e9; % Carrier Frequency fc = 2.4 GHz
sample_freq_before = 20e6; % Sampling Rate (before upsampling) = 20 MHz
sample_freq_after = 5e9; % Sampling Rate (after upsampling) = 5 GHz


%% Part 1 OFDM Signal Generation and Reception
bps = 4; % 4 bits per symbol
M = 2^bps; % Modulation Order: 16-QAM

txsymbols = randi([0 M - 1], fft_size, 1);
txmod = qammod(txsymbols, M, UnitAveragePower=true); % Use 16-QAM
txout = ifft(txmod, fft_size);

txout = txout(:);
txcp = txout(fft_size - cyclic_prefix_size + 1 : fft_size);
txout = [txcp; txout];

% 1A. Generate OFDM symbol
hold on
stem(txout, 'filled', 'b');
stem(txcp, 'filled', 'r')
title("OFDM symbol with cyclic prefix")
xlabel("time (10^{-6} s)")
legend("output signal", "cyclic prefix")

% 1B. OFDM symbol after upconversion in frequency domain


% 1C. Plot received signal in frequency domain

% 1D. Transmitting 1000 OFDM symbols with different SNR

% 1E. Repeat 1D using 64-QAM
%% Part 2 STS & LTS Generation and Frame Synchronization

% 2A. Generate STS
STS_unit = sqrt(13/6)*[0, 0, 1 + 1j, 0, 0, 0, -1 - 1j, 0, 0, 0, 1 + 1j, 0, 0, 0, -1 - 1j ,0, 0, 0, -1 - 1j, 0, 0, 0, 1 + 1j, 0, 0, 0, 0, 0, 0, 0, -1 - 1j, 0, 0, 0, -1 - 1j, 0, 0, 0, 1 + 1j, 0, 0, 0, 1 + 1j, 0, 0, 0, 1 + 1j, 0, 0, 0, 1 + 1j, 0, 0];
STS = zeros(1, 10*53);

for i = 1:10
    STS( (53*(i - 1) + 1):(53*i) ) = STS_unit;
end

figure;
stem(STS, 'filled', 'b')
title("STS sequence")
xlabel("Time (10^{-6} s)")


% 2B. Generate LTS
LTS_unit = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
LTS = zeros(1, 2 * 53);

for i = 1:2
    LTS( (53 * (i - 1) + 1): (53 * i) ) = LTS_unit;
end

LTS = LTS(:);
LTS_cp = LTS(81:106);
LTS = [LTS_cp; LTS];

figure;
hold on
stem(LTS, 'filled','b')
stem(LTS_cp, 'filled', 'r')
title("LTS sequence")
xlabel("Time (10^{-6} s)")
legend("LTS", "Cyclic Prefix")

% 2C. Frame with 3 OFDM signals

% 2D. Frame with 100 OFDM signals

% 2E. Explain Synchronization Algorithms

% 2F. BER of signal in 2D

% 2G. Transmission with Different SNR
    
