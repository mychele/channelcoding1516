%% Viterbi encoder & decoder simulator
rng default

%% Param
N = 3;
R = 1/N;
rho = 2/3;
pck_len = 10^4;
it_num = 10^4;
ebn0_dB = -3:8;
ebn0 = 10.^(ebn0_dB/10);
sigma_w = sqrt(1./(2/3*ebn0));
mem = 2; % code memory
max_err = 10^3;

%% Simulation of Viterbi, SD
n_err_sd_ci = zeros(size(ebn0));
pck_sent_sd_ci = zeros(size(ebn0));

tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	u = randi(2, pck_len + mem, 1) - 1;
	u(pck_len + 1:pck_len + mem) = 0;
	% encode
	c = encoder517_mex(u.');
	% conform map
	s = 2*c.' - 1;
	% random gaussian noise
	w = randn(length(c), 1);
	for j = 1:length(sigma_w)
		if n_err_sd_ci(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w(j), 0);
			n_err_sd_ci(j) = n_err_sd_ci(j) + sum(u ~= u_hat.');
			pck_sent_sd_ci(j) = pck_sent_sd_ci(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_ci./(pck_sent_sd_ci*pck_len))]);
	end
end
time_sd_ci = toc
save('517_BER.mat');

%% Simulation of Viterbi, SD, simpler input
n_err_sd_si = zeros(size(ebn0));
pck_sent_sd_si = zeros(size(ebn0));

tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w)
		if n_err_sd_si(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w(j), 0);
			n_err_sd_si(j) = n_err_sd_si(j) + sum(u_hat.' ~= 0);
			pck_sent_sd_si(j) = pck_sent_sd_si(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_si./(pck_sent_sd_si*pck_len))]);
	end
end
time_sd_si = toc
save('517_BER.mat');

%% Simulation of Viterbi, SD, simpler input, windowed 3
n_err_sd_win3 = zeros(size(ebn0));
pck_sent_sd_win3 = zeros(size(ebn0));

winSize = 3;
tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w)
		if n_err_sd_win3(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex_win(y.', sigma_w(j), 0, winSize);
			n_err_sd_win3(j) = n_err_sd_win3(j) + sum(u_hat.' ~= 0);
			pck_sent_sd_win3(j) = pck_sent_sd_win3(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win3./(pck_sent_sd_win3*pck_len))]);
	end
end
time_sd_win3 = toc
save('517_BER.mat');


%% Simulation of Viterbi, SD, simpler input, windowed 5
n_err_sd_win5 = zeros(size(ebn0));
pck_sent_sd_win5 = zeros(size(ebn0));

winSize = 5;
tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w)
		if n_err_sd_win5(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex_win(y.', sigma_w(j), 0, winSize);
			n_err_sd_win5(j) = n_err_sd_win5(j) + sum(u_hat.' ~= 0);
			pck_sent_sd_win5(j) = pck_sent_sd_win5(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win5./(pck_sent_sd_win5*pck_len))]);
	end
end
time_sd_win5 = toc
save('517_BER.mat');

%% Simulation of Viterbi, SD, simpler input, windowed 10
n_err_sd_win10 = zeros(size(ebn0));
pck_sent_sd_win10 = zeros(size(ebn0));

winSize = 10;
tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w)
		if n_err_sd_win10(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex_win(y.', sigma_w(j), 0, winSize);
			n_err_sd_win10(j) = n_err_sd_win10(j) + sum(u_hat.' ~= 0);
			pck_sent_sd_win10(j) = pck_sent_sd_win10(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win10./(pck_sent_sd_win10*pck_len))]);
	end
end
time_sd_win10 = toc
save('517_BER.mat');


%% Simulation of Viterbi, HD, simpler input
n_err_hd = zeros(size(ebn0));
pck_sent_hd = zeros(size(ebn0));

tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w)
		if n_err_hd(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w(j), 1);
			n_err_hd(j) = n_err_hd(j) + sum(u_hat.' ~= 0);
			pck_sent_hd(j) = pck_sent_hd(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_hd./(pck_sent_hd*pck_len))]);
	end
end
time_hd = toc
save('517_BER.mat');

%% Compute nominal and effective coding gain
ebn0_dB = -3:13;
Pbit_coded = @(x) qfunc(sqrt(4*x)) + 4*qfunc(sqrt(16/3*x)) + 12*qfunc(sqrt(20/3*x));
Pbit_uncoded = @(x) qfunc(sqrt(2*x));
gamma_c = R*6;
Pbit_nom = qfunc(sqrt(2*10.^(ebn0_dB/10)*gamma_c));

Pe_target = 10^-5;
diff_target_uncoded = @(x) abs(Pbit_uncoded(x) - Pe_target);
ebno_uncoded = fminsearch(@(x) diff_target_uncoded(x), 10);
diff_target_coded = @(x) abs(Pbit_coded(x) - Pe_target);
ebno_coded = fminsearch(@(x) diff_target_coded(x), 10);
disp(['Nominal coding gain = ' num2str(10*log10(gamma_c)) ' dB'])
disp(['Effective coding gain = ' num2str(10*log10(ebno_uncoded/ebno_coded)) ' dB'])
disp(['Difference = ' num2str(10*log10(gamma_c) - 10*log10(ebno_uncoded/ebno_coded)) ' dB'])

figure,
semilogy(ebn0_dB, Pbit_coded(10.^(ebn0_dB/10)), '-d', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_nom, '-s', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_uncoded(10.^(ebn0_dB/10)), '-o', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('EB/N0')
ylabel('BER')
ylim([10^-10, 1])
grid on
legend('Bound ()', 'Expression ()', 'Uncoded')

%% Plot simulation against bounds
linewidth = 1.2;
load('517_BER.mat');

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*pck_len), '-x', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_si./(pck_sent_sd_si*pck_len), '-o', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_coded(10.^(ebn0_dB/10)), '-d', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_nom, '-s', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_uncoded(10.^(ebn0_dB/10)), '-o', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('$\frac{E_b}{N_0}$')
ylabel('BER')
grid on
legend('Simulation, random input', 'Simulation, zero input', '\eqref{eq:BER_bound}', '\eqref{eq:BER_approx}', 'Uncoded')



% semilogy(ebn0_dB, n_err_sd_win3./(pck_sent_sd_win3*pck_len), '-^'), hold on
% semilogy(ebn0_dB, n_err_sd_win5./(pck_sent_sd_win5*pck_len), '-v'), hold on
% semilogy(ebn0_dB, n_err_sd_win10./(pck_sent_sd_win10*pck_len), '-*'), hold on
% semilogy(ebn0_dB, n_err_hd./(pck_sent_hd*pck_len), '-.'), hold on






