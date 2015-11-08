%% a plot (use Matlab) of the expected performance (Pbit versus Eb/N0)
% together with the derivation of both the nominal and the effective coding gain

ebn0 = 10.^((0:10)./10);
R = 1/3;
Pbit_code = qfunc(sqrt(4*ebn0)) + 4*qfunc(sqrt(16/3*ebn0)) + 12*qfunc(sqrt(20/3*ebn0));
Pbit_uncoded = qfunc(sqrt(2*ebn0));

gamma_c = R*6;
Pbit_nom = qfunc(sqrt(2*ebn0*gamma_c));

%load('517_BER.mat');

BER = sum(BER_mat, 2)/10;

figure, 
plot(10*log10(ebn0), log10(Pbit_code)), hold on
plot(10*log10(ebn0), log10(Pbit_uncoded)), hold on
plot(10*log10(ebn0), log10(Pbit_nom)), hold on
plot(ebn0_dB, log10(BER)), hold on
legend('g3', 'uncoded', 'nominal coding gain', 'viterbi')

