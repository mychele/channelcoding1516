%% Viterbi encoder & decoder simulator
rng default

%% Param
N = 3;
rho = 2/3;
pck_len = 10^4;
it_num = 10^4;
ebn0_dB = 0:8;
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

%% Plot
Pbit_code = qfunc(sqrt(4*ebn0)) + 4*qfunc(sqrt(16/3*ebn0)) + 12*qfunc(sqrt(20/3*ebn0));
Pbit_uncoded = qfunc(sqrt(2*ebn0));
gamma_c = R*6;
Pbit_nom = qfunc(sqrt(2*ebn0*gamma_c));

figure,
plot(ebn0_dB, log10(n_err_sd_ci./(pck_sent_sd_ci*pck_len)), '-x'), hold on
plot(ebn0_dB, log10(n_err_sd_si./(pck_sent_sd_si*pck_len)), '-o'), hold on
plot(ebn0_dB, log10(n_err_sd_win3./(pck_sent_sd_win3*pck_len)), '-^'), hold on
plot(ebn0_dB, log10(n_err_sd_win5./(pck_sent_sd_win5*pck_len)), '-v'), hold on
plot(ebn0_dB, log10(n_err_sd_win10./(pck_sent_sd_win10*pck_len)), '-*'), hold on
plot(ebn0_dB, log10(n_err_hd./(pck_sent_hd*pck_len)), '-.'), hold on
plot(ebn0_dB, log10(Pbit_code), '--'), hold on
plot(ebn0_dB, log10(Pbit_uncoded), ':'), hold on
plot(ebn0_dB, log10(Pbit_nom), '-s'), hold on
xlabel('EB/N0')
ylabel('BER')
legend('sd, full input', 'sd, 0 input', 'sd, win3', 'sd, win5', 'sd, win10', 'hd')






