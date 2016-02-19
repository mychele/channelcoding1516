close all

K = 30592;
N = 32640;

ebno_dB = [4.5, 5, 5.5, 6, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9];
gamma = 10.^(ebno_dB/10)*K/N;
gamma_dB = 10*log10(gamma);

BER = [0.0107698, 0.00746274, 0.00497123, 0.00316586, 0.00197241, 0.00145332, 0.000851121, 0.000366599, 7.93508e-05, 4.3312e-06];

Pbit_uncoded = @(x) qfunc(sqrt(2*x));

figure,
semilogy(gamma_dB, BER, '-^'), hold on,
semilogy(0:10, Pbit_uncoded(10.^((0:10)/10)), '-v')
xlim([0, 10])
grid on
xlabel('Gamma [dB]')
ylabel('BER')
legend('LDPC, simulation', 'Uncoded')

