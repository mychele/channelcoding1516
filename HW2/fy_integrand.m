function [ integrand ] = fy_integrand( b, sigma_w, a, px )
% This function computes the integrand of the differential entropy for 
% M-ary PAM modulation with soft decoding

fy_vec = 1/sqrt(2*pi*sigma_w^2) * exp(-(b-a).^2/(2*sigma_w^2)).'*px;

integrand = -sum(fy_vec)*log2(sum(fy_vec));

if (isnan(integrand) || isinf(integrand))
        integrand = 0;
end 

end

