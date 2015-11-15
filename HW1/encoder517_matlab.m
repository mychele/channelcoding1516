function [ c ] = encoder517_matlab( u )
% This function encodes a random input u into a 517 coded output
% @param the input vector, column
% @return the codeword vector, column

N = 3; % number of filters
mem = 2; % memory of the code
state = [0; 0];
c = zeros(length(u)*N, 1);

v1 = [1; 0; 1];
v2 = [1; 0];
m1 = [0 1; 0 1; 1 1];
m2 = [0 0; 1 0];

for l = 0:length(u)-1
	c(l*N+1:(l+1)*N) = mod(v1*u(l+1) + m1*state, 2);
	state = mod(v2*u(l+1) + m2*state, 2);
end

