%% a plot (use Matlab) of the expected performance (Pbit versus Eb/N0)
% together with the derivation of both the nominal and the effective coding gain

d = [6; 8; 10];
Kd = [1; 4; 12];
snr = 10.^((0:15)./10);
R = 1/3;
Pbit_d = bsxfun(@times, Kd, qfunc(sqrt(2*R*d*snr)));
Pbit_code = sum(Pbit_d, 1);
Pbit_uncoded = qfunc(sqrt(2*snr));

gamma_c = R*min(d);
Pbit_nok = qfunc(sqrt(2*snr*gamma_c));

open('517_BER.mat');

BER = mean(BER_mat, 2);

figure, 
plot(10*log10(snr), log10(Pbit_code)), hold on
plot(10*log10(snr), log10(Pbit_uncoded)), hold on
plot(10*log10(snr), log10(Pbit_nok)), hold on
plot(10*log10(snr), log10(BER)), hold on
legend('g3', 'uncoded', 'nominal coding gain', 'viterbi')

