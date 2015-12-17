close all;
clear;
clc

%% Data for plot
color_matrix = colormap;
markers = {'+-','o-','*-','x-','s-','d-','^-','v-','>-','<-','p-','h-'};
markersize = 4;
linewidth = 1.2;

%% Integration for mutual information - uniform input distribution
snr_dB = 0:40;
M_vec = [2,4,8,16,32,64];
inforate_uni = zeros(length(snr_dB), length(M_vec));

for M_ind = 1:length(M_vec)
    % params
    M = M_vec(M_ind);
    disp(M);
    
    px = 1/M*ones(M, 1); % column vector that describes the px(a)
    d = 0:M-1;
    d = d(:);
    a = -(M-1) + 2*d;
    a = a(:); % make it a column
    
    % compute avg symbol energy
    es = px.' * abs(a).^2;
    
    for snr_in = 1:length(snr_dB)
        sigma_w = sqrt( es /  10^( snr_dB(snr_in)/10 )  );
        alpha = -0.5*log2(2*pi*exp(1)*sigma_w^2);
        integrand = @(b) fy_integrand(b, sigma_w, a, px);
        hy = integral(integrand, -Inf, Inf, 'ArrayValued', 1);
        inforate_uni(snr_in, M_ind) = hy + alpha;
    end
end

figure, hold on
for i = 1:length(M_vec)
    line_fewer_markers(snr_dB, 2*inforate_uni(:, i), 21, markers{mod(i, numel(markers)) + 1}, ...
        'Color', color_matrix(mod(i*10, size(color_matrix, 1)) + 1,:), ...
        'LineWidth', linewidth, 'DisplayName', strcat('M= ', num2str(M_vec(i))), 'MarkerSize', markersize)
end
i = i +1;
line_fewer_markers(snr_dB, log2(1+10.^(snr_dB/10)), 21, markers{mod(i, numel(markers)) + 1}, ...
        'Color', color_matrix(mod(i*10, size(color_matrix, 1)) + 1,:), ...
        'LineWidth', linewidth, 'DisplayName', 'Shannon bound', 'MarkerSize', markersize)
legend('-DynamicLegend', 'Location', 'NorthWest')
grid on
xlabel('SNR [dB]')
ylabel('2I [bit/s/Hz]')

%% Gaussian input distribution and maximization

snr_dB = 0:40;
M_vec = [2,4,8,16,32,64];
inforate_gauss = zeros(length(snr_dB), length(M_vec));
sigma_in_max = zeros(length(snr_dB), length(M_vec));

for snr_in = 1:length(snr_dB)
    
    disp(snr_dB(snr_in))
    for M_ind = 1:length(M_vec)
        % params
        M = M_vec(M_ind);
        d = 0:M-1;
        d = d(:);
        a = -(M-1) + 2*d;
        a = a(:); % make it a column
        sigma_in_vec = (0.4:0.2:10)*M;
        
        inforate_snr_vec = zeros(length(sigma_in_vec), 1);
        
        for sigma_in_ind = 1:length(sigma_in_vec)
            sigma_in = sigma_in_vec(sigma_in_ind);
            px = 1/(sqrt(2*pi)*sigma_in)*exp(-a.^2/(2*sigma_in^2)); % column vector that describes the px(a)
            % normalize px
            px = px/sum(px);
            % compute avg symbol energy
            es = px.' * abs(a).^2;  
            
            sigma_w = sqrt( es /  10^( snr_dB(snr_in)/10 )  );
            alpha = -0.5*log2(2*pi*exp(1)*sigma_w^2);
            
            integrand = @(b) fy_integrand(b, sigma_w, a, px);
            hy = integral(integrand, -Inf, Inf, 'ArrayValued', 1);
            inforate_snr_vec(sigma_in_ind) = hy + alpha;
        end
        [inforate_gauss(snr_in, M_ind), sigma_max] = max(inforate_snr_vec);
        sigma_in_max(snr_in, M_ind) = sigma_in_vec(sigma_max);
    end
    
end

figure, hold on
for i = 1:length(M_vec)
    line_fewer_markers(snr_dB, 2*inforate_gauss(:, i), 21, markers{mod(i, numel(markers)) + 1}, ...
        'Color', color_matrix(mod(i*12, size(color_matrix, 1)) + 1,:), ...
        'LineWidth', linewidth, 'DisplayName', strcat('M= ', num2str(M_vec(i))), 'MarkerSize', markersize)
end
i = i+1;
line_fewer_markers(snr_dB, log2(1+10.^(snr_dB/10)), 21, markers{mod(i, numel(markers)) + 1}, ...
        'Color', color_matrix(mod(i*12, size(color_matrix, 1)) + 1,:), ...
        'LineWidth', linewidth, 'DisplayName', 'Shannon bound', 'MarkerSize', markersize)
legend('-DynamicLegend', 'Location', 'NorthWest')
grid on
xlabel('SNR [dB]')
ylabel('2C {bit/s/Hz]')


%% More refined search for the maximum sigma_IN for each SNR
snr_dB = 20:35;
M = 32;
inforate_gauss = zeros(length(snr_dB), length(M_vec));
sigma_in_max = zeros(length(snr_dB), 1);

for snr_in = 1:length(snr_dB)
    disp(snr_dB(snr_in))
    % params
    d = 0:M-1;
    d = d(:);
    a = -(M-1) + 2*d;
    a = a(:); % make it a column
    sigma_in_vec = 0:160;
    
    inforate_snr_vec = zeros(length(sigma_in_vec), 1);
    
    for sigma_in_ind = 1:length(sigma_in_vec)
        sigma_in = sigma_in_vec(sigma_in_ind);
        px = 1/(sqrt(2*pi)*sigma_in)*exp(-a.^2/(2*sigma_in^2)); % column vector that describes the px(a)
        
        % normalize px
        px = px/sum(px);
        
        % compute avg symbol energy
        es = px.' * abs(a).^2;
        sigma_w = sqrt( es /  10^( snr_dB(snr_in)/10 )  );
        alpha = -0.5*log2(2*pi*exp(1)*sigma_w^2);
        
        integrand = @(b) fy_integrand(b, sigma_w, a, px);
        hy = integral(integrand, -Inf, Inf, 'ArrayValued', 1);
        inforate_snr_vec(sigma_in_ind) = hy + alpha;
        
    end
    
    [inforate_gauss(snr_in), sigma_max] = max(inforate_snr_vec);
    sigma_in_max(snr_in) = sigma_in_vec(sigma_max);
    
end

figure,
color_matrix = colormap;
markers = {'+-','o-','*-','x-','s-','d-','^-','v-','>-','<-','p-','h-'};
markersize = 4;
linewidth = 1.2;
i = 1;
line_fewer_markers(snr_dB, sigma_in_max, 21, markers{mod(i, numel(markers)) + 1}, ...
        'Color', color_matrix(mod(i*12, size(color_matrix, 1)) + 1,:), ...
        'LineWidth', linewidth, 'DisplayName', '\sigma_{IN}', 'MarkerSize', markersize)
legend('-DynamicLegend', 'Location', 'NorthWest')
grid on
xlabel('SNR [dB]')
ylabel('\sigma_{IN}')



