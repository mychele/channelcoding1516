sigma_w2 = 1;
xmax = 10;
x = -xmax:xmax; % support (received signal)
x_plus = 0:xmax;
x_min = -xmax:0;
R = 1;
Pbit = qfunc(sqrt(1/sigma_w2*R));

%% AWGN, soft decoding
f_awgnsd = 2*x/sigma_w2;

%% AWGN, hard decoding
f_awgnhd_plus = log((1-Pbit)/Pbit)*ones(1, length(x_plus));
f_awgnhd_minus = -log((1-Pbit)/Pbit)*ones(1, length(x_min)-1);
f_awgnhd = [f_awgnhd_minus, f_awgnhd_plus];


%% Laplacian
f_lapl = -sqrt(2)/sqrt(sigma_w2)*(abs(x-1) - abs(x+1));

%% MoG
p_x1 = sqrt(5)/(8*sqrt(2*pi*sigma_w2))*(4/sqrt(2)*exp(-5*(x-1).^2/(4*sigma_w2)) + exp(-5*(x-1).^2/(8*sigma_w2)) + ...
	 1/sqrt(8)*exp(-5*(x-1).^2/(16*sigma_w2)) + 1/4*exp(-5*(x-1).^2/(32*sigma_w2)));
p_x0 = sqrt(5)/(8*sqrt(2*pi*sigma_w2))*(4/sqrt(2)*exp(-5*(x+1).^2/(4*sigma_w2)) + exp(-5*(x+1).^2/(8*sigma_w2)) + ...
	 1/sqrt(8)*exp(-5*(x+1).^2/(16*sigma_w2)) + 1/4*exp(-5*(x+1).^2/(32*sigma_w2)));
f_mog = log(p_x1./p_x0);

%% plot
figure,
plot(x, f_awgnsd, x, f_lapl, x, f_mog), hold on
stairs(x, f_awgnhd)
xlabel('x')
ylabel('LLR')
legend('AWGN, SD', 'LAPLACE', 'MoG', 'AWGN, HD')