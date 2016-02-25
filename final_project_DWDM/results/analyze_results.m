close all

K = 30592;
N = 32640;

Pbit_uncoded = @(x) qfunc(sqrt(2*x));
PER_uncoded = @(x) 1 - (1 - Pbit_uncoded(x)).^K;

%% Read simulation_results_multi files for SumProduct
ebno_dB_sumproduct = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.25, 6.5, 6.6, 6.7:0.05:7.1];
gamma_sumproduct = 10.^(ebno_dB_sumproduct/10)*K/N;
gamma_dB_sumproduct = 10*log10(gamma_sumproduct);
num_err_packet_sp = zeros(length(ebno_dB_sumproduct), 1);
num_err_bit_sp = zeros(length(ebno_dB_sumproduct), 1);
num_attempt_per_snr_sp = zeros(length(ebno_dB_sumproduct), 1);

for i = 1:length(ebno_dB_sumproduct)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/simulation_results_multi_', num2str(ebno_dB_sumproduct(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    num_attempt_per_snr_sp(i) = length(vec);
    num_err_packet_sp(i) = length(vec(vec~=0));
    num_err_bit_sp(i) = sum(vec(vec~=0));
end

BER_sp = num_err_bit_sp./num_attempt_per_snr_sp;
BER_sp = BER_sp./K;

PER_sp = num_err_packet_sp./num_attempt_per_snr_sp;

ebno_dB_minsum = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.25, 6.5, 6.6, 6.7:0.05:7.15];
gamma_minsum = 10.^(ebno_dB_minsum/10)*K/N;
gamma_dB_minsum = 10*log10(gamma_minsum);
num_err_packet_ms = zeros(length(ebno_dB_minsum), 1);
num_err_bit_ms = zeros(length(ebno_dB_minsum), 1);
num_attempt_per_snr_ms = zeros(length(ebno_dB_minsum), 1);

for i = 1:length(ebno_dB_minsum)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/simulation_results_minsum', num2str(ebno_dB_minsum(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    num_attempt_per_snr_ms(i) = length(vec);
    num_err_packet_ms(i) = length(vec(vec~=0));
    num_err_bit_ms(i) = sum(vec(vec~=0));
end

BER_ms = num_err_bit_ms./num_attempt_per_snr_ms;
BER_ms = BER_ms./K;

PER_ms = num_err_packet_ms./num_attempt_per_snr_ms;

figure,
semilogy(gamma_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(gamma_dB_minsum, BER_ms, '-*'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([4, 10])
ylim([5e-8, 3e-1])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('LDPC, simulation, SumProduct', 'LDPC, simulation, MinSum', 'Uncoded')


figure,
semilogy(gamma_dB_sumproduct, PER_sp, '-^'), hold on,
semilogy(gamma_dB_minsum, PER_ms, '-*'), hold on,
semilogy(0:10, PER_uncoded(10.^((0:10)/10)), '-v')
xlim([4, 10])
ylim([3e-5, 1.05])
grid on
xlabel('Gamma [dB]')
ylabel('PER')
legend('LDPC, simulation, SumProduct', 'LDPC, simulation, MinSum', 'Uncoded')


%% MinSum

time_minsum_n = [1.68831e+09, 1.24743e+09, 7.90332e+08, 5.26343e+08, 3.4851e+08, 1.98934e+08, 1.49833e+08, 1.35503e+08, 1.24461e+08, 1.17984e+08];
time_minsum = time_minsum_n./(1e9); % in second

time_sumproduct_n = [1.03241e+10, 9.12977e+09, 7.99721e+09, 7.10846e+09, 6.55603e+09, 4.30983e+09, 3.80277e+09, 3.52308e+09, 3.3068e+09];
time_sumproduct = time_sumproduct_n./(1e9); % in second


%% Time
figure,
stem(gamma_dB_sumproduct(end-8:end), time_sumproduct, '-^'), hold on,
stem(gamma_dB_minsum(end-9:end), time_minsum, '-*'), hold on,
grid on
xlabel('Gamma [dB]')
ylabel('Time to perform a decoding [s]')
legend('LDPC, SumProduct', 'LDPC, MinSum')



