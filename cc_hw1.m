%% Viterbi encoder & decoder simulator
rng default

%% Param
N = 3;
rho = 2/3;
pck_len = 10^4;
it_num = 10^5;
ebn0_dB = 0:10;
ebn0 = 10.^(ebn0_dB/10);
sigma_w = sqrt(1./(2/3*ebn0));
mem = 2; % code memory
max_err = 10^3;
n_err = zeros(size(ebn0));
pck_sent = zeros(size(ebn0));

%% Simulation
tic
for it_index = 1:it_num
	% define input vector of size pck_len + mem
	u = randi(2, pck_len + mem, 1) - 1;
	u(pck_len + 1:pck_len + mem) = 0;
	% encode
	c = encoder517_matlab(u);
	% conform map
	s = 2*c - 1;
	% random gaussian noise
	w = randn(length(c), 1);
	
	for j = 1:length(sigma_w)
		if n_err(j) < max_err
			% channel output
			y = s + sigma_w(j)*w;
			% decode with Viterbi
			u_hat = viterbi_mex(y.', sigma_w(j));
			n_err(j) = n_err(j) + sum(u ~= u_hat.');
			pck_sent(j) = pck_sent(j) + 1;
		end
	end	
	if(mod(it_index, 1000) == 0)
		save('517_BER_mex.mat');
		% display current status
		toc
        disp(['#' num2str(it_index) ', BER = ' num2str(n_err./(pck_sent*pck_len))]);
	end
end