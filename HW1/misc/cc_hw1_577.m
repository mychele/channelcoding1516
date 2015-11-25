%% Param
N = 3;
R = 1/N;
rho = 2/3;
pck_len = 10^4;
it_num = 5*10^5;
ebn0_dB = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8];
ebn0 = 10.^(ebn0_dB/10);
sigma_w = sqrt(1./(2/3*ebn0));
mem = 2; % code memory
max_err = 10^3;
save_data = 1;
save_int = 10000;

%% Simulation of Viterbi, SD
n_err_sd_ci_577 = zeros(size(ebn0));
pck_sent_sd_ci_577 = zeros(size(ebn0));

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

%% Plot 517 vs 577
linewidth = 1.2;

load('517_BER.mat');
load('577_BER.mat');

figure,
semilogy(ebn0_dB, n_err_sd_ci./(pck_sent_sd_ci*pck_len), '-x', 'LineWidth', linewidth, 'MarkerSize', 4), hold on
semilogy(ebn0_dB, n_err_sd_ci_577./(pck_sent_sd_ci_577*pck_len), '-o', 'LineWidth', linewidth, 'MarkerSize', 4), hold on


xlabel('$\frac{E_b}{N_0}$')
ylabel('BER')
grid on
legend('517, random input', '577, random input')

