N = 3; % codeword length
mem = 2; % memory length (size of the state vector)
M = 2; % binary constellation
num_states = M^mem;
mu = 10;

% Precompute all the possible output vectors
y_lut = zeros(M*N, num_states); % each column i has 2 codeword y, from state 
% i with input 0 or 1. Each codeword has size N
state = zeros(mem, 1); % state vector of size mem, with lsb in the left most position!
state_ID = 0:3; % state ID (i.e. 00 -> 0, 01 -> 1...)
input_vec_mult = [1; 0; 1];
state_mat_mult = [0 1; 0 1; 1 1];
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

% Neighbors vector
neigh = zeros(num_states, mem);
for state_ID = 0:num_states-1
	switch mod(state_ID, 2)
		case 0
			neigh(state_ID + 1, :) = [0, 1];
		case 1
			neigh(state_ID + 1, :) = [2, 3];
			
	end
end

% Cycle on the received vector r
% Init
Gamma_prev = inf*ones(num_states, 1);
Gamma_prev(1) = 0; % start from state 0
Gamma = inf*ones(num_states, 1);

prev_state = zeros(mu+mem, num_states);

while (Gamma(1) == inf)
	% cycle through states
	for state_ID = 0:num_states-1
		% compute the input symbol that can bring to this state
		% for 00 and 01 the input symbol that bring to the state is always
		% 0, for 10 and 11 the input symbol is 1
		u_poss = floor(state_ID/mem);
		
		% cycle on the neighbors, compute new cost, update if lower than
		% before
		cost = inf;
		for neigh_ID = neigh(state_ID+1, :)
			newcost = Gamma_prev(neigh_ID+1) + sum(y_lut(u_poss*N + 1:(u_poss + 1)*N, neigh_ID+1));
			if (newcost < cost) % new min
				if(~(state_ID==0 && neigh_ID == 0 && newcost == 0)) % avoid to account for the path from 0 to 0
					cost = newcost;
					prev_state(ceil(l/N), state_ID+1) = neigh_ID;
				end
				Gamma(state_ID+1) = cost;
			end	
		end
		if (cost == inf) % it was not initialized
			prev_state(ceil(l/N), state_ID+1) = neigh(1);
			Gamma(state_ID+1) = cost;
		end
	end
	Gamma_prev = Gamma;
	l = l + 1;
end

