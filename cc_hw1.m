%% Viterbi encoding and decoding
%rng default
% decode with mex -- in C

N = 3;
rho = 2/3;
%Gamma_dB = 1:8;
%Gamma = 10.^(Gamma_dB/10);
pck_len = 10^4;
it_num = 10;
ebn0_dB = 1:8;
ebn0 = 10.^(ebn0_dB/10);
sigma_w = sqrt(1./(2/3*ebn0));

mem = 2; % code memory
BER_mat = zeros(length(ebn0), it_num);
% 
% u = randi(2, pck_len + mem, 1) - 1;
% u(pck_len + 1:pck_len + mem) = 0;
% 
% % encode
% c = encoder517_matlab(u);
% 
% % conform map
% s = 2*c - 1;
% 
% % random gaussian noise
% w = randn(length(c), 1);
% 
% % channel output
% y = s + sigma_w(1)*w;
% 
% % decode with Viterbi
% u_hat = viterbi517_matlab(y, sigma_w(1), u, s);
% 
% figure,
% stem(u), hold on,
% stem(u_hat),
% legend('u', 'u_hat')

for i = 1:it_num
	disp(i);
	tic;
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
		% channel output
		y = s + sigma_w(j)*w;
		
		% decode with Viterbi
		u_hat = viterbi517_matlab(y, sigma_w(j));
		BER_mat(j, i) = sum(u ~= u_hat)/pck_len;
	end
	toc
	
	if(mod(i, 200) == 0)
		save('517_BER.mat');
	end
end