%% a plot (use Matlab) of the expected performance (Pbit versus Eb/N0)
% together with the derivation of both the nominal and the effective coding gain

load('517_BER_mex.mat');
BER_mex = n_err./(pck_sent*pck_len);
prec = 1 - log10(pck_len*it_num);
disp(['Expected precision = 10^' num2str(prec)])

load('517_BER_matlab.mat');
BER_matlab = mean(BER_mat, 2);
prec = 1 - log10(pck_len*it_num);
disp(['Expected precision = 10^' num2str(prec)])

ebn0 = 10.^((0:9)./10);
R = 1/3;
Pbit_code = qfunc(sqrt(4*ebn0)) + 4*qfunc(sqrt(16/3*ebn0)) + 12*qfunc(sqrt(20/3*ebn0));
Pbit_uncoded = qfunc(sqrt(2*ebn0));

gamma_c = R*6;
Pbit_nom = qfunc(sqrt(2*ebn0*gamma_c));


figure, 
plot(10*log10(ebn0), log10(Pbit_code), '.-'), hold on
plot(10*log10(ebn0), log10(Pbit_uncoded)), hold on
plot(10*log10(ebn0), log10(Pbit_nom)), hold on
plot(ebn0_dB(1:end-1), log10(BER_mex(1:end-1)), '^-'), hold on
plot(ebn0_dB(1:end-2), log10(BER_matlab(1:end-2)), 'v-'), hold on

legend('g3', 'uncoded', 'nominal coding gain', 'viterbi mex', 'viterbi matlab')

