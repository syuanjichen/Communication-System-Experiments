%% Question 1. Compare Natural Code v.s. Gray Code Using 16-QAM

M_QAM16 = 16;      % Modulation order
k_QAM16 = log2(M_QAM16); % Number of bits per symbol
n = 4500000;  % Number of symbols per frame
sps = 1;     % Number of samples per symbol (oversampling factor)
dB_points = 80; % Number of points of Eb/No
rng default  % Use default random number generator

dataIn = randi([0 1], n * k_QAM16, 1); % Generate vector of binary data

dataSymbolsIn = bit2int(dataIn, k_QAM16);

dataMod_QAM16_bin = qammod(dataSymbolsIn, M_QAM16, 'bin'); % Binary-encoded
dataMod_QAM16 = qammod(dataSymbolsIn, M_QAM16);      % Gray-encoded

decibel = linspace(0, 20, dB_points + 1);
ber_EbNo_QAM16_bin = zeros(1, dB_points + 1);
ber_EbNo_QAM16 = zeros(1, dB_points + 1);
ber_EbNo_BPSK = zeros(1, dB_points + 1);
ber_EbNo_QPSK = zeros(1, dB_points + 1);
ber_EbNo_PSK8 = zeros(1, dB_points + 1);
ber_EbNo_PSK16 = zeros(1, dB_points + 1);

for i = 1 : dB_points + 1
    EbNo = 10^(0.1 * decibel(i));
    snr = convertSNR(EbNo,'ebno', ...
        samplespersymbol = sps, ...
        bitspersymbol = k_QAM16);
    
    receivedSignal_QAM16_bin = awgn(dataMod_QAM16_bin, snr, 'measured');
    receivedSignal_QAM16 = awgn(dataMod_QAM16, snr, 'measured');
    
    dataSymbolsOut_QAM16_bin = qamdemod(receivedSignal_QAM16_bin , M_QAM16, 'bin'); % Binary-encoded data symbols
    dataSymbolsOut_QAM16 = qamdemod(receivedSignal_QAM16, M_QAM16);     % Gray-coded data symbols
    
    dataOut_QAM16_bin = int2bit(dataSymbolsOut_QAM16_bin, k_QAM16);
    dataOut_QAM16 = int2bit(dataSymbolsOut_QAM16, k_QAM16);
    
    [numErrors_QAM16_bin, ber_16QAM_bin] = biterr(dataIn, dataOut_QAM16_bin);
    ber_EbNo_QAM16_bin(i) = ber_16QAM_bin;
    fprintf('\n i = %d, The binary coding bit error rate is %5.2e, based on %d errors.\n', ...
       i, ber_16QAM_bin, numErrors_QAM16_bin)
    
    [numErrors_QAM16, ber_QAM16] = biterr(dataIn,dataOut_QAM16);
    ber_EbNo_QAM16(i) = ber_QAM16;
    fprintf('\n i = %d, The Gray coding bit error rate is %5.2e, based on %d errors.\n', ...
        i, ber_QAM16, numErrors_QAM16)
end

ber_EbNo_QAM16_theo = berawgn(decibel, 'qam', M_QAM16);

figure
plot(decibel, ber_EbNo_QAM16_theo, '--o', 'LineWidth', 2)
hold on
plot(decibel, ber_EbNo_QAM16_bin, '-xr', 'LineWidth', 2)
hold on
plot(decibel, ber_EbNo_QAM16, '-.b', 'LineWidth', 2)
hold off
title('16-QAM: Binary Code vs Gray Code')
xlabel('E_b / N_0 (dB)')
ylabel('Bit Error Rate')
yscale log
legend('Theo', 'Real: Natural Code', 'Real: Gray Code')

%% Question 2. Compare BPSK, QPSK, and 8-PSK Using Gray Code

M_BPSK = 2;
k_BPSK = log2(M_BPSK);

M_QPSK = 4;
k_QPSK = log2(M_QPSK);

M_PSK8 = 8;
k_PSK8 = log2(M_PSK8);

dataSymbolsIn_BPSK = bit2int(dataIn, k_BPSK);
dataSymbolsIn_QPSK = bit2int(dataIn, k_QPSK);
dataSymbolsIn_PSK8 = bit2int(dataIn, k_PSK8);

dataMod_BPSK = pskmod(dataSymbolsIn_BPSK, M_BPSK);
dataMod_QPSK = pskmod(dataSymbolsIn_QPSK, M_QPSK);
dataMod_PSK8 = pskmod(dataSymbolsIn_PSK8, M_PSK8);

for i = 1 : dB_points + 1
    EbNo = 10^(0.1 * decibel(i));

    snr_BPSK = convertSNR(EbNo,'ebno', ...
        samplespersymbol = sps, ...
        bitspersymbol = 1);

    snr_QPSK = convertSNR(EbNo,'ebno', ...
        samplespersymbol = sps, ...
        bitspersymbol = 2);

    snr_8PSK = convertSNR(EbNo,'ebno', ...
        samplespersymbol = sps, ...
        bitspersymbol = 3);
    
    receivedSignal_BPSK = awgn(dataMod_BPSK, snr_BPSK, 'measured');
    receivedSignal_QPSK = awgn(dataMod_QPSK, snr_QPSK, 'measured');
    receivedSignal_PSK8 = awgn(dataMod_PSK8, snr_8PSK, 'measured');

    dataSymbolsOut_BPSK = pskdemod(receivedSignal_BPSK, M_BPSK);
    dataSymbolsOut_QPSK = pskdemod(receivedSignal_QPSK, M_QPSK);
    dataSymbolsOut_PSK8 = pskdemod(receivedSignal_PSK8, M_PSK8);
    
    dataOut_BPSK = int2bit(dataSymbolsOut_BPSK, 1);
    dataOut_QPSK = int2bit(dataSymbolsOut_QPSK, 2);
    dataOut_PSK8 = int2bit(dataSymbolsOut_PSK8, 3);
    
    [numErrors_BPSK, ber_BPSK] = biterr(dataIn, dataOut_BPSK);
    ber_EbNo_BPSK(i) = ber_BPSK;
    fprintf('\n i = %d, The BPSK bit error rate is %5.2e, based on %d errors.\n', ...
       i, ber_BPSK, numErrors_BPSK)
    
    [numErrors_QPSK, ber_QPSK] = biterr(dataIn, dataOut_QPSK);
    ber_EbNo_QPSK(i) = ber_QPSK;
    fprintf('\n i = %d, The QPSK bit error rate is %5.2e, based on %d errors.\n', ...
        i, ber_QPSK, numErrors_QPSK)

    [numErrors_PSK8, ber_PSK8] = biterr(dataIn, dataOut_PSK8);
    ber_EbNo_PSK8(i) = ber_PSK8;
    fprintf('\n i = %d, The 8-PSK bit error rate is %5.2e, based on %d errors.\n', ...
        i, ber_PSK8, numErrors_PSK8)
end

ber_EbNo_theo_BPSK = berawgn(decibel, 'psk', 2, 'nondiff');
ber_EbNo_theo_QPSK = berawgn(decibel, 'psk', 4, 'nondiff');
ber_EbNo_theo_PSK8 = berawgn(decibel, 'psk', 8, 'nondiff');

figure;
plot(decibel, ber_EbNo_theo_BPSK, '--^', 'LineWidth', 2)
hold on
plot(decibel, ber_EbNo_BPSK, '-xr', 'LineWidth', 2)
plot(decibel, ber_EbNo_theo_QPSK, '-*', 'LineWidth', 2)
plot(decibel, ber_EbNo_QPSK, '-ob', 'LineWidth', 2)
plot(decibel, ber_EbNo_theo_PSK8, '--o', 'LineWidth', 2)
plot(decibel, ber_EbNo_PSK8, '-square', 'LineWidth', 2)
hold off
yscale log
title('BPSK vs QPSK vs 8-PSK')
xlabel('E_b / N_0 (dB)')
ylabel('Bit Error Rate')
legend('Theo BPSK', 'Real BPSK', 'Theo QPSK', 'Real QPSK', 'Theo 8-PSK', 'Real 8-PSK')

%% Question 3. Compare 16-PSK and 16-QAM Using Gray Code

M_PSK16 = 16;
k_PSK16 = log2(M_PSK16);

dataSymbolsIn_PSK16 = bit2int(dataIn, k_PSK16);

dataMod_PSK16 = pskmod(dataSymbolsIn_PSK16, M_PSK16);

for i = 1 : dB_points + 1
    EbNo = 10^(0.1 * decibel(i));

    snr_16bit = convertSNR(EbNo,'ebno', ...
        samplespersymbol = sps, ...
        bitspersymbol = 4);
    
    receivedSignal_PSK16 = awgn(dataMod_PSK16, snr_16bit, 'measured');

    dataSymbolsOut_PSK16 = pskdemod(receivedSignal_PSK16, 16);
    
    dataOut_PSK16 = int2bit(dataSymbolsOut_PSK16, 4);
    
    [numErrors_PSK16, ber_PSK16] = biterr(dataIn, dataOut_PSK16);
    ber_EbNo_PSK16(i) = ber_PSK16;
    fprintf('\n i = %d, The 16-PSK bit error rate is %5.2e, based on %d errors.\n', ...
       i, ber_PSK16, numErrors_PSK16)
end

ber_EbNo_PSK16_theo = berawgn(decibel, 'psk', 16, 'nondiff');

figure;
hold on
plot(decibel, ber_EbNo_PSK16_theo, '--o', 'LineWidth', 2)
plot(decibel, ber_EbNo_PSK16, '-xr', 'LineWidth', 2)
plot(decibel, ber_EbNo_QAM16_theo, '--o', 'LineWidth', 2)
plot(decibel, ber_EbNo_QAM16, '-b', 'marker', '.', 'LineWidth', 2)
hold off

yscale log
title('16-PSK vs 16-QAM')
xlabel('E_b / N_0 (dB)')
ylabel('Bit Error Rate')
legend('Theo 16-PSK', 'Real 16-PSK', 'Theo 16-QAM', 'Real 16-QAM')
