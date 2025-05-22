clear all;
clc;

%% Part 1 Successful Transmission with Up/Down Sampling

num_bits = 10; % Number of bits
modulation_order = 2; % BPSK is used, so order = 2
bits_per_symbol = log2(modulation_order); % A symbol contains log2(M) bits
L = 6; % Interpolation Factor
M2 = 3; % Decimation Factor
rng default % Use Default Random Number Generator

omega = linspace(-2*pi, 2*pi, num_bits); % Frequency samples

data_bits_in = randi([0, 1], 1, num_bits); % Create a random input with size num_bits

%stem(data_bits_in) % Plot input signal in Time-Domain

data_mod_BPSK = pskmod(data_bits_in, modulation_order); % Modulate the Signals

data_mod_upsampled = upsample(data_mod_BPSK, L); % Upsampling the Modulated Signal

b1 = fir1(48, 1/L);
%data_mod_interp = filter(b1, 1, data_mod_upsampled);
data_mod_interp = lowpass(data_mod_upsampled, 1/L);
%plot(omega_upsample, abs(fftshift(fft(data_mod_interp))));

b2 = fir1(48, 1/M2);
%data_mod_decimate = filter(b2, 1, data_mod_interp);
data_mod_decimate = lowpass(data_mod_interp, 1/M2);
%plot(omega_upsample, abs(fftshift(fft(data_mod_decimate))));

%data_mod_decimate = decimate(data_mod_intfilt, L);
data_mod_downsampled = downsample(data_mod_decimate, M2);

data_demod = pskdemod(data_mod_downsampled, modulation_order);
data_bits_out = downsample(data_demod, L/M2);

%stem(data_bits_out)

[numErrors, ber] = biterr(data_bits_in, data_bits_out);
fprintf("Case 1, Input Size: %d, Errors: %d, Error rate = %.3f \n", num_bits, numErrors, ber)

% Plot the time domain signals
figure;
subplot(2, 1, 1)
hold on
stem(1:num_bits, data_mod_BPSK, 'filled', ':diamondb')
stem(1:1/L:num_bits + 1 - 1/L, data_mod_upsampled, 'filled')
hold off
xlabel("n")
legend("BPSK", "Upsampled")
subplot(2, 1, 2)
hold on
stem(1:1/L:num_bits + 1 - 1/L, data_mod_interp, 'filled', ':diamondb')
stem(1:M2/L:num_bits + 1 - M2/L, data_mod_downsampled, 'filled')
hold off
xlabel("n")
legend("Interpolated", "Decimated")

% Plot the frequency domain signals
figure;
hold on
plot(-2*pi:4*pi/length(data_mod_BPSK):2*pi - 4*pi/length(data_mod_BPSK), abs(fft(fftshift(data_mod_BPSK))))
plot(-2*pi:4*pi/length(data_mod_upsampled):2*pi - 4*pi/length(data_mod_upsampled), abs(fft(fftshift(data_mod_upsampled))))
plot(-2*pi:4*pi/length(data_mod_interp):2*pi - 4*pi/length(data_mod_interp), abs(fft(fftshift(data_mod_interp))))
hold off
xlabel("\omega (rad/s)")
ylabel("|X(e^{j\omega)})|")
legend("BPSK", "Upsampled", "Interpolated")

%% Part 2 Aliasing in Downsampling with or without the Decimator Filter

data_bits_in_2 = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0];
M2 = 10;
num_bits = length(data_bits_in_2);

data_mod_BPSK_2 = pskmod(data_bits_in_2, modulation_order); % Modulate the Signals
data_mod_upsampled_2 = upsample(data_mod_BPSK_2, L); % Upsampling the Modulated Signal

b1 = fir1(48, 1/L);
%data_mod_interp = filter(b1, 1, data_mod_upsampled);
data_mod_interp_2 = lowpass(data_mod_upsampled_2, 1/L);

b2 = fir1(48, 1/M2);
%data_mod_decimate = filter(b2, 1, data_mod_interp);
data_mod_decimate_2 = lowpass(data_mod_interp_2, 1/M2);
%plot(omega_upsample, abs(fftshift(fft(data_mod_decimate))));

%data_mod_decimate = decimate(data_mod_intfilt, L);
data_mod_downsampled_nodecimate_2 = downsample(data_mod_interp_2, M2);
data_mod_downsampled_2 = downsample(data_mod_decimate_2, M2);

data_demod_2 = pskdemod(data_mod_downsampled_2, modulation_order);
data_demod_nodecimate_2 = pskdemod(data_mod_downsampled_nodecimate_2, modulation_order);
%data_bits_out_2 = downsample(data_demod_2, L/M2);
%data_bits_out_nodecimate_2 = downsample(data_demod_nodecimate_2, L/M2);

%stem(data_bits_out)

[numErrors_2, ber_2] = biterr(data_bits_in_2, data_bits_out_2);
[numErrors_nodecimate_2, ber_nodecimate_2] = biterr(data_bits_in_2, data_bits_out_nodecimate_2);
fprintf("Case 2, With Decimation, Input Size: %d, Errors: %d, Error rate = %.3f \n", num_bits, numErrors_2, ber_2)
fprintf("Case 2, Without Decimation, Input Size: %d, Errors: %d, Error rate = %.3f \n", num_bits, numErrors_nodecimate_2, ber_nodecimate_2)

% Plot the time domain signals
figure;
subplot(2, 1, 1)
hold on
stem(1:num_bits, data_mod_BPSK_2, 'filled', ':diamondb')
stem(1:1/L:num_bits + 1 - 1/L, data_mod_upsampled_2, 'filled')
hold off
xlabel("n")
legend("BPSK", "Upsampled")
subplot(2, 1, 2)
hold on
stem(1:1/L:num_bits + 1 - 1/L, data_mod_interp_2, 'filled', ':diamondb')
stem(1:M2/L:num_bits + 1 - M2/L, data_mod_downsampled_2, 'filled')
stem(1:M2/L:num_bits + 1 - M2/L, data_mod_downsampled_nodecimate_2, 'filled')
hold off
xlabel("n")
legend("Interpolated", "No Decimation", "Decimated")

% Plot the frequency domain signals
figure;
hold on
plot(-2*pi:4*pi/length(data_mod_BPSK_2):2*pi - 4*pi/length(data_mod_BPSK_2), abs(fft(fftshift(data_mod_BPSK_2))))
plot(-2*pi:4*pi/length(data_mod_upsampled_2):2*pi - 4*pi/length(data_mod_upsampled_2), abs(fft(fftshift(data_mod_upsampled_2))))
plot(-2*pi:4*pi/length(data_mod_interp_2):2*pi - 4*pi/length(data_mod_interp_2), abs(fft(fftshift(data_mod_interp_2))))
hold off
xlabel("\omega (rad/s)")
ylabel("|X(e^{j\omega)})|")
legend("BPSK", "Upsampled", "Interpolated")

% Change Input
data_bits_in_3 = [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0];

data_mod_BPSK_3 = pskmod(data_bits_in_3, modulation_order); % Modulate the Signals
data_mod_upsampled_3 = upsample(data_mod_BPSK_3, L); % Upsampling the Modulated Signal

b1 = fir1(48, 1/L);
%data_mod_interp = filter(b1, 1, data_mod_upsampled);
data_mod_interp_3 = lowpass(data_mod_upsampled_3, 1/L);

b2 = fir1(48, 1/M2);
%data_mod_decimate = filter(b2, 1, data_mod_interp);
data_mod_decimate_3 = lowpass(data_mod_interp_3, 1/M2);
%plot(omega_upsample, abs(fftshift(fft(data_mod_decimate))));

%data_mod_decimate = decimate(data_mod_intfilt, L);
data_mod_downsampled_nodecimate_3 = downsample(data_mod_interp_3, M2);
data_mod_downsampled_3 = downsample(data_mod_decimate_3, M2);

data_demod_3 = pskdemod(data_mod_downsampled_3, modulation_order);
data_demod_nodecimate_3 = pskdemod(data_mod_downsampled_nodecimate_3, modulation_order);
data_bits_out_3 = downsample(data_demod_3, L/M2);
data_bits_out_nodecimate = downsample(data_demod_nodecimate_3, L/M2);

% Plot the time domain signals
figure;
subplot(2, 1, 1)
hold on
stem(1:num_bits, data_mod_BPSK_3, 'filled', ':diamondb')
stem(1:1/L:num_bits + 1 - 1/L, data_mod_upsampled_3, 'filled')
hold off
xlabel("n")
legend("BPSK", "Upsampled")
subplot(2, 1, 2)
hold on
stem(1:1/L:num_bits + 1 - 1/L, data_mod_interp_3, 'filled', ':diamondb')
stem(1:M2/L:num_bits + 1 - M2/L, data_mod_downsampled_3, 'filled')
stem(1:M2/L:num_bits + 1 - M2/L, data_mod_downsampled_nodecimate_3, 'filled')
hold off
xlabel("n")
legend("Interpolated", "No Decimation", "Decimated")

%% Part 3 Square Root Raised Cosine (RRC) Filter 

alpha = 0.5; % Roll-off Factor of Square Root Raised Cosine Filter
    
