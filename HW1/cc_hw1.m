%% Viterbi encoder & decoder simulator
rng default

%% Param
N = 3;
R = 1/N;
rho = 2/3;
pck_len = 10^4;
it_num = 10^4;
ebn0_dB = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8];
ebn0 = 10.^(ebn0_dB/10);
sigma_w = sqrt(1./(2/3*ebn0));
mem = 2; % code memory
max_err = 10^3;
save_data = 0;
save_int = 1000;

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
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat');
		end
		% display current status
		time_sd_ci(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_ci./(pck_sent_sd_ci*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
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
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat');
		end
		% display current status
		time_sd_si(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_si./(pck_sent_sd_si*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
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
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat');
		end		% display current status
		time_win3(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win3./(pck_sent_sd_win3*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
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
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat');
		end		% display current status
		time_win5(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win5./(pck_sent_sd_win5*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
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
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat');
		end		% display current status
		time_win10(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_win10./(pck_sent_sd_win10*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
%% Simulation of Viterbi, HD, simpler input
ebn0_dB_hdsi = -3:12;
n_err_hd_si = zeros(size(ebn0_dB_hdsi));
pck_sent_hd_si = zeros(size(ebn0_dB_hdsi));
sigma_w_hdsi = sqrt(1./(2/3*(10.^(ebn0_dB_hdsi/10))));

tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	s = -1*ones(N*(pck_len+mem), 1);
	% random gaussian noise
	w = randn(length(s), 1);
	
	for j = 1:length(sigma_w_hdsi)
		if n_err_hd_si(j) < max_err
			% channel output
			y = s + sigma_w_hdsi(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w_hdsi(j), 1);
			n_err_hd_si(j) = n_err_hd_si(j) + sum(u_hat.' ~= 0);
			pck_sent_hd_si(j) = pck_sent_hd_si(j) + 1;
		end
	end	
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('517_BER.mat', '-append');
		end		% display current status
		time_hd(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_hd_si./(pck_sent_hd_si*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
%% Simulation of Viterbi, HD, complete input
ebn0_dB_hdci = -3:8;
n_err_hd_ci = zeros(size(ebn0_dB_hdci));
pck_sent_hd_ci = zeros(size(ebn0_dB_hdci));
sigma_w_hdci = sqrt(1./(2/3*(10.^(ebn0_dB_hdci/10))));

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
	
	for j = 1:length(sigma_w_hdci)
		if n_err_hd_ci(j) < max_err
			% channel output
			y = s + sigma_w_hdci(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w_hdci(j), 1);
			n_err_hd_ci(j) = n_err_hd_ci(j) + sum(u ~= u_hat.');
			pck_sent_hd_ci(j) = pck_sent_hd_ci(j) + 1;
		end
	end	
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			
			save('517_BER.mat', '-append');
		end		% display current status
		time_hd_ci(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_hd_ci./(pck_sent_hd_ci*(pck_len+mem)))]);
	end
end

if(save_data == 1)
	save('517_BER.mat');
end
%% Simulation of Viterbi, SD, 577
n_err_sd_ci_577 = zeros(size(ebn0));
pck_sent_sd_ci_577 = zeros(size(ebn0));
%it_num = 5*10^5;

tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	u = randi(2, pck_len + mem, 1) - 1;
	u(pck_len + 1:pck_len + mem) = 0;
	% encode
	c = encoder577_mex(u.');
	% conform map
	s = 2*c.' - 1;
	% random gaussian noise
	w = randn(length(c), 1);
	for j = 1:length(sigma_w)
		if n_err_sd_ci_577(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex577(y.', sigma_w(j), 0);
			n_err_sd_ci_577(j) = n_err_sd_ci_577(j) + sum(u ~= u_hat.');
			pck_sent_sd_ci_577(j) = pck_sent_sd_ci_577(j) + 1;
		end
	end	
	if(mod(it_index, save_int) == 0)
		if(save_data == 1)
			save('577_BER.mat');
		end
		% display current status
		time_sd_ci_577(it_index/save_int) = toc;
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err_sd_ci_577./(pck_sent_sd_ci_577*pck_len))]);
	end
end

if(save_data == 1)
	save('577_BER.mat');
end

%% Elaborate time vectors
load('ex_time_100.mat');
int_sd_ci(2:length(time_sd_ci)) = time_sd_ci(2:end) - time_sd_ci(1:end-1);
int_sd_ci(1) = time_sd_ci(1);
int_sd_si(2:length(time_sd_si)) = time_sd_si(2:end) - time_sd_si(1:end-1);
int_sd_si(1) = time_sd_si(1);
int_win3(2:length(time_win3)) = time_win3(2:end) - time_win3(1:end-1);
int_win3(1) = time_win3(1);
int_win5(2:length(time_win5)) = time_win5(2:end) - time_win5(1:end-1);
int_win5(1) = time_win5(1);
int_win10(2:length(time_win10)) = time_win10(2:end) - time_win10(1:end-1);
int_win10(1) = time_win10(1);
int_hd(2:length(time_hd)) = time_hd(2:end) - time_hd(1:end-1);
int_hd(1) = time_hd(1);
int_hd_ci(2:length(time_hd_ci)) = time_hd_ci(2:end) - time_hd_ci(1:end-1);
int_hd_ci(1) = time_hd_ci(1);

mean_sd_ci = mean(int_sd_ci);
ci_sd_ci = 1.96*std(int_sd_ci)/sqrt(length(int_sd_ci));
mean_sd_si = mean(int_sd_si);
ci_sd_si = 1.96*std(int_sd_si)/sqrt(length(int_sd_si));
mean_win3 = mean(int_win3);
ci_win3 = 1.96*std(int_win3)/sqrt(length(int_win3));
mean_win5 = mean(int_win5);
ci_win5 = 1.96*std(int_win5)/sqrt(length(int_win5));
mean_win10 = mean(int_win10);
ci_win10 = 1.96*std(int_win10)/sqrt(length(int_win10));
mean_hd = mean(int_hd);
ci_hd = 1.96*std(int_hd)/sqrt(length(int_hd));
mean_hd_ci = mean(int_hd_ci);
ci_hd_ci = 1.96*std(int_hd_ci)/sqrt(length(int_hd_ci));

mean_vec = [mean_sd_ci, mean_sd_si, mean_win3, mean_win5, mean_win10, mean_hd, mean_hd_ci];
ci_vec = [ci_sd_ci, ci_sd_si, ci_win3, ci_win5, ci_win10, ci_hd, ci_hd_ci];
figure,
errorbar(mean_vec, ci_vec, 'o'),
grid on,
xlabel('Different simulation setup')
ylabel('Time to execute 100 simulations [s]')

%% Compute nominal and effective coding gain
linewidth = 1.2;
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
semilogy(ebn0_dB, Pbit_coded(10.^(ebn0_dB/10)), '-d', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_nom, '-s', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_uncoded(10.^(ebn0_dB/10)), '-o', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('EB/N0')
ylabel('BER')
ylim([10^-10, 1])
grid on
legend('Bound ()', 'Expression ()', 'Uncoded')

%% Plot simulation against bounds
load('517_BER.mat');
Pbit_nom = qfunc(sqrt(2*10.^(ebn0_dB/10)*gamma_c));

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*(pck_len+mem)), '-x', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_si./(pck_sent_sd_si*(pck_len+mem)), '-o', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_coded(10.^(ebn0_dB/10)), '-d', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_nom, '-s', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, Pbit_uncoded(10.^(ebn0_dB/10)), '-o', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('$\frac{E_b}{N_0}$')
ylabel('BER')
grid on
legend('Simulation, random input', 'Simulation, zero input', ...
	'\eqref{eq:BER_bound}', '\eqref{eq:BER_approx}', 'Uncoded')

%% Plot simulation against windowed
load('517_BER.mat');
Pbit_nom = qfunc(sqrt(2*10.^(ebn0_dB/10)*gamma_c));

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*(pck_len+mem)), '-x', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_si./(pck_sent_sd_si*(pck_len+mem)), '-o',...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_win3./(pck_sent_sd_win3*(pck_len+mem)), '-^',...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_win5./(pck_sent_sd_win5*(pck_len+mem)), '-v', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_win10./(pck_sent_sd_win10*(pck_len+mem)), '-*', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('\frac{E_b}{N_0}')
ylabel('BER')
grid on
legend('Simulation, random input', 'Simulation, zero input', 'Windowed, zero input, \upsilon = 3', ...
	'Windowed, zero input, \upsilon = 5', 'Windowed, zero input, \upsilon = 10')

%% Plot simulation against HD
load('517_BER.mat');
Pbit_nom = qfunc(sqrt(2*10.^(ebn0_dB/10)*gamma_c));

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*(pck_len+mem)), '-x', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_si./(pck_sent_sd_si*(pck_len+mem)), '-o', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB_hdsi, n_err_hd_si./(pck_sent_hd_si*(pck_len+mem)), ...
	'-*', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB_hdci, n_err_hd_ci./(pck_sent_hd_ci*(pck_len+mem)), ...
	'-^', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB_hdci, Pbit_uncoded(10.^(ebn0_dB_hdci/10)), '--', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('\frac{E_b}{N_0}')
ylabel('BER')
xlim([-3, 10])
ylim([10^-7, 1])
grid on
legend('SD, random input', 'SD, zero input', 'HD, zero input', ...
	'HD, random input', 'Uncoded')

%% Plot 517 vs 577
load('517_BER.mat');
load('577_BER.mat');

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*pck_len), '-x', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_ci_577./(pck_sent_sd_ci_577*pck_len), '-o', ...
	'LineWidth', linewidth, 'MarkerSize', 4), hold on
xlabel('$\frac{E_b}{N_0}$')
ylabel('BER')
grid on
legend('517, random input', '577, random input')