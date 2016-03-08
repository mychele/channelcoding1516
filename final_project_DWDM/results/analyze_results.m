close all

K = 30592;
N = 32640;

Pbit_uncoded = @(x) qfunc(sqrt(2*x));
PER_uncoded = @(x) 1 - (1 - Pbit_uncoded(x)).^K;

%% Read simulation_results_multi files for SumProduct
ebno_dB_sumproduct = [3, 3.5, 4, 4.2, 4.6, 4.7:0.02:4.94];
gamma_sumproduct = 10.^(ebno_dB_sumproduct/10)*K/N*2;
gamma_dB_sumproduct = 10*log10(gamma_sumproduct);
num_err_packet_sp = zeros(length(ebno_dB_sumproduct), 1);
num_err_bit_sp = zeros(length(ebno_dB_sumproduct), 1);
num_attempt_per_snr_sp = zeros(length(ebno_dB_sumproduct), 1);

for i = 1:length(ebno_dB_sumproduct)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/blade50/simulation_results_multi_', num2str(ebno_dB_sumproduct(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_sp(i) = length(vec);
    num_err_packet_sp(i) = length(vec(vec~=0));
    num_err_bit_sp(i) = sum(vec(vec~=0));
end

BER_sp = num_err_bit_sp./num_attempt_per_snr_sp;
BER_sp = BER_sp./K;
PER_sp_fromBER = 1 - (1 - BER_sp).^K;
PER_sp = num_err_packet_sp./num_attempt_per_snr_sp;

figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([2, 10])
ylim([5e-8, 1e-1])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('LDPC', 'Uncoded')

figure,
semilogy(ebno_dB_sumproduct, PER_sp, '-^'), hold on,
xlim([2, 7])
ylim([5e-6, 1.1])
grid on
xlabel('EBNO [dB]')
ylabel('PER')
legend('50 iterations')

%% Read simulation results for 1 iteration
ebno_dB_1sp = [2,3,4,4.5,5,5.5,6,6.5,7];
gamma_1sp = 10.^(ebno_dB_1sp/10)*K/N;
gamma_dB_1sp = 10*log10(gamma_1sp);
num_err_packet_1sp = zeros(length(ebno_dB_1sp), 1);
num_err_bit_1sp = zeros(length(ebno_dB_1sp), 1);
num_attempt_per_snr_1sp = zeros(length(ebno_dB_1sp), 1);

for i = 1:length(ebno_dB_1sp)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/blade1/simulation_results_multi_1_', num2str(ebno_dB_1sp(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_1sp(i) = length(vec);
    num_err_packet_1sp(i) = length(vec(vec~=0));
    num_err_bit_1sp(i) = sum(vec(vec~=0));
end

BER_1sp = num_err_bit_1sp./num_attempt_per_snr_1sp;
BER_1sp = BER_1sp./K;
PER_1sp_fromBER = 1 - (1 - BER_1sp).^K;

figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_1sp, BER_1sp, '-*'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([2, 7])
ylim([1e-8, 1e-1])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('50', '1', 'Uncoded')

%% Read results for 10 iterations
ebno_dB_10sp = [2,3,4,4.5,4.6,4.7,4.8,4.9,5,5.1,6];
gamma_10sp = 10.^(ebno_dB_10sp/10)*K/N;
gamma_dB_10sp = 10*log10(gamma_10sp);
num_err_packet_10sp = zeros(length(ebno_dB_10sp), 1);
num_err_bit_10sp = zeros(length(ebno_dB_10sp), 1);
num_attempt_per_snr_10sp = zeros(length(ebno_dB_10sp), 1);

for i = 1:length(ebno_dB_10sp)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/blade10/simulation_results_multi_10_', num2str(ebno_dB_10sp(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_10sp(i) = length(vec);
    num_err_packet_10sp(i) = length(vec(vec~=0));
    num_err_bit_10sp(i) = sum(vec(vec~=0));
end

BER_10sp = num_err_bit_10sp./num_attempt_per_snr_10sp;
BER_10sp = BER_10sp./K;
PER_10sp_fromBER = 1 - (1 - BER_10sp).^K;

figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_10sp, BER_10sp, '-*'), hold on,
semilogy(ebno_dB_1sp, BER_1sp, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([2, 7])
ylim([1e-8, 1e-1])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('50', '10', '1', 'Uncoded')

%% Read results for 20 iterations
ebno_dB_20sp = [3,4,4.2,4.5,4.6,4.7,4.8,4.9, 4.94, 5];
gamma_20sp = 10.^(ebno_dB_20sp/10)*K/N;
gamma_dB_20sp = 10*log10(gamma_20sp);
num_err_packet_20sp = zeros(length(ebno_dB_20sp), 1);
num_err_bit_20sp = zeros(length(ebno_dB_20sp), 1);
num_attempt_per_snr_20sp = zeros(length(ebno_dB_20sp), 1);

for i = 1:length(ebno_dB_20sp)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/blade20/simulation_results_multi_20_', num2str(ebno_dB_20sp(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_20sp(i) = length(vec);
    num_err_packet_20sp(i) = length(vec(vec~=0));
    num_err_bit_20sp(i) = sum(vec(vec~=0));
end

BER_20sp = num_err_bit_20sp./num_attempt_per_snr_20sp;
BER_20sp = BER_20sp./K;
PER_20sp_fromBER = 1 - (1 - BER_20sp).^K;

figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_20sp, BER_20sp, '-d'), hold on,
semilogy(ebno_dB_10sp, BER_10sp, '-*'), hold on,
semilogy(ebno_dB_1sp, BER_1sp, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([4.3, 8])
ylim([1e-8, 1e-1])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('50', '20', '10', '1', 'Uncoded')

%% Results for 100 iterations
ebno_dB_100sp = [4.5,4.6,4.65,4.7:0.02:4.94];
gamma_100sp = 10.^(ebno_dB_100sp/10)*K/N;
gamma_dB_100sp = 10*log10(gamma_100sp);
num_err_packet_100sp = zeros(length(ebno_dB_100sp), 1);
num_err_bit_100sp = zeros(length(ebno_dB_100sp), 1);
num_attempt_per_snr_100sp = zeros(length(ebno_dB_100sp), 1);

for i = 1:length(ebno_dB_100sp)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/blade100/simulation_results_multi_100_', num2str(ebno_dB_100sp(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_100sp(i) = length(vec);
    num_err_packet_100sp(i) = length(vec(vec~=0));
    num_err_bit_100sp(i) = sum(vec(vec~=0));
end

BER_100sp = num_err_bit_100sp./num_attempt_per_snr_100sp;
BER_100sp = BER_100sp./K;
PER_100sp_fromBER = 1 - (1 - BER_100sp).^K;

figure,
semilogy(ebno_dB_100sp, BER_100sp, '-x'), hold on,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_20sp, BER_20sp, '-d'), hold on,
semilogy(ebno_dB_10sp, BER_10sp, '-*'), hold on,
semilogy(ebno_dB_1sp, BER_1sp, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([4.5, 4.95])
ylim([5e-8, 2e-2])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('100', '50', '20', '10', '1', 'Uncoded')

figure,
semilogy(ebno_dB_100sp, BER_100sp, '-x'), hold on,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_20sp, BER_20sp, '-d'), hold on,
semilogy(ebno_dB_10sp, BER_10sp, '-*'), hold on,
semilogy(ebno_dB_1sp, BER_1sp, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([3.5, 6])
ylim([5e-8, 2e-2])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('100', '50', '20', '10', '1', 'Uncoded')

%% Compare with MinSum


ebno_dB_minsum =  [2, 2.5, 3, 3.5, 4, 4.5, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5]; 
gamma_minsum = 10.^(ebno_dB_minsum/10)*K/N;
gamma_dB_minsum = 10*log10(gamma_minsum);
num_err_packet_ms = zeros(length(ebno_dB_minsum), 1);
num_err_bit_ms = zeros(length(ebno_dB_minsum), 1);
num_attempt_per_snr_ms = zeros(length(ebno_dB_minsum), 1);

for i = 1:length(ebno_dB_minsum)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/bladems50/simulation_results_minsum_', num2str(ebno_dB_minsum(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_ms(i) = length(vec);
    num_err_packet_ms(i) = length(vec(vec~=0));
    num_err_bit_ms(i) = sum(vec(vec~=0));
end

BER_ms = num_err_bit_ms./num_attempt_per_snr_ms;
BER_ms = BER_ms./K;
PER_ms_fromBER = 1 - (1 - BER_ms).^K;
PER_ms = num_err_packet_ms./num_attempt_per_snr_ms;


figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_minsum, BER_ms, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([4, 6])
ylim([5e-8, 2e-2])
grid on
xlabel('EBN0 [dB]')
ylabel('BER')
legend('Sum Product', 'Min Sum', 'Uncoded')


%% Different slopes
ebno_dB_ossp = [2,3,4,4.5,4.6,4.7,4.74,4.78,4.8:0.02:4.94];
gamma_ossp = 10.^(ebno_dB_ossp/10)*K/N;
gamma_dB_ossp = 10*log10(gamma_ossp);
num_err_packet_ossp = zeros(length(ebno_dB_ossp), 1);
num_err_bit_ossp = zeros(length(ebno_dB_ossp), 1);
num_attempt_per_snr_ossp = zeros(length(ebno_dB_ossp), 1);

for i = 1:length(ebno_dB_ossp)
    fileID = fopen(strcat('/Users/Michele/Dev/git/channelcoding1516/final_project_DWDM/results/results_blade/bladeos/simulation_results_os100_', num2str(ebno_dB_ossp(i)), '.txt'), 'r');
    vec = fscanf(fileID, '%d');
    fclose(fileID);
    num_attempt_per_snr_ossp(i) = length(vec);
    num_err_packet_ossp(i) = length(vec(vec~=0));
    num_err_bit_ossp(i) = sum(vec(vec~=0));
end

BER_ossp = num_err_bit_ossp./num_attempt_per_snr_ossp;
BER_ossp = BER_ossp./K;
PER_ossp_fromBER = 1 - (1 - BER_ossp).^K;

figure,
semilogy(ebno_dB_sumproduct, BER_sp, '-^'), hold on,
semilogy(ebno_dB_ossp, BER_ossp, '-o'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([4, 5.5])
ylim([5e-8, 2e-2])
grid on
xlabel('EBN0 [dB]')
ylabel('BER')
legend('1,2,3,4,5,6,7', '1,2,3,5,7,11,13', 'Uncoded')

% 
ebno_sp = [3, 4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5];
time_sp_50 = [6.21475e+09, 6.17265e+09, 6.03251e+09, 5.04901e+09, 2.60614e+09,  1.39746e+09, 1.07678e+09, 6.03046e+08, 3.72513e+08, 3.41108e+08, 3.05817e+08, 2.82775e+08, 2.66328e+08]./1e9;
time_sp_100 = [1.06808e+10, 1.05435e+10, 1.01902e+10, 7.64292e+09, 2.64309e+09, 1.26275e+09, 9.9787e+08, 6.39012e+08, 3.78606e+08, 3.42002e+08, 3.14643e+08, 2.85593e+08, 2.67953e+08]./1e9;

ebno_ms = [3, 4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5];
time_ms = [4.18477e+09, 4.19191e+09, 4.1868e+09, 4.20782e+09,  4.18253e+09, 4.20487e+09, 4.19902e+09, 4.01526e+09, 3.42146e+09, 1.77495e+09, 7.43431e+08, 4.72207e+08, 3.99272e+08]./1e9;


figure,
stem(1:length(time_sp_50), log10(time_sp_50), '-^'), hold on,
stem(1:length(time_sp_100), log10(time_sp_100), '-x'), hold on,
stem(1:length(time_ms),  log10(time_ms), '-*'), hold on,
grid on
xlabel('EBN0 [dB]')
ylabel('T [s]')
legend('Sum Product, 50', 'Sum Product, 100', 'Min Sum, 50')

