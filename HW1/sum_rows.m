function [ H ] = sum_rows( H, src, dst )
% Sum src to dst in H (operations in F2) and returns H
H(dst, :) = mod(H(src, :) + H(dst, :), 2);
end

