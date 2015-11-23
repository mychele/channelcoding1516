N = 2; % codeword length
mem = 2; % memory length (size of the state vector)
M = 2; % binary constellation
num_states = M^mem;
mu = 20;

% Precompute all the possible output vectors
y_lut = zeros(M*N, num_states); % each column i has 2 codeword y, from state 
% i with input 0 or 1. Each codeword has size N
state = zeros(mem, 1); % state vector of size mem, with lsb in the left most position!
state_ID = 0:3; % state ID (i.e. 00 -> 0, 01 -> 1...)

% g3
input_vec_mult = [1; 0; 1];
state_mat_mult = [0 1; 0 1; 1 1];

% g4
% input_vec_mult = [1; 0; 1];
% state_mat_mult = [0 0; 0 1; 1 0];

% (5, 7)
% input_vec_mult = [1; 1];
% state_mat_mult = [0 1; 1 1];

% (5, 7, 7)
% input_vec_mult = [1; 1; 1];
% state_mat_mult = [0 1; 1 1; 1 1];

% state_ID, symbol
%neigh_ext = [0, 0; 1, 1; 2, 0; 3, 1; 0, 1; 1, 0; 2, 1; 3, 0]; % g4
%neigh_ext = [0, 0; 1, 0; 2, 0; 3, 0; 0, 1; 1, 1; 2, 1; 3, 1]; % g3 but also (5,7)
neigh_ext = [0, 0; 1, 1; 2, 1; 3, 0; 0, 1; 1, 0; 2, 0; 3, 1];

input_vec_mult = [1; 1];
state_mat_mult = [0 0; 1 0];

for state_ID = 0:num_states-1
	% update state vector
	state(1) = mod(state_ID, mem);
	state(2) = floor(state_ID/mem);
	
	out_state_part = state_mat_mult*flipud(state);
	for input_symbol = 0:1
		out_in_part = input_vec_mult.*input_symbol;
		y_lut(input_symbol*N + 1:(input_symbol + 1)*N, state_ID+1) = mod(out_in_part + out_state_part, 2);
	end
end

% Cycle on the received vector r
% Init
Gamma_prev = inf*ones(num_states, 1);
Gamma_prev(1) = 0; % start from state 0
Gamma = inf*ones(num_states, 1);

for l = 1:20
	% cycle through states
	for state_ID = 0:num_states-1
		% cycle on the neighbors, compute new cost, update if lower than
		% before
		cost = inf;
		for i = 1:mem
			neigh_ID = neigh_ext(state_ID*mem + i, 1);
			u_poss = neigh_ext(state_ID*mem + i, 2);
			if (neigh_ID == 0 && u_poss == 0)
				newcost = inf; % invalidate transition from 0 to 0
			else
				newcost = Gamma_prev(neigh_ID+1) + sum(y_lut(u_poss*N + 1:(u_poss + 1)*N, neigh_ID+1));
			end
			if (newcost < cost) % new min
				if(~(state_ID==0 && neigh_ID == 0 && newcost == 0)) % avoid to account for the path from 0 to 0
					cost = newcost;
				end
				Gamma(state_ID+1) = cost;
			end
		end
		if (cost == inf) % it was not initialized
			Gamma(state_ID+1) = cost;
		end
	end
	Gamma_prev = Gamma;
 	Gamma_zero(l) = Gamma(1);
end
disp(min(Gamma_zero))
