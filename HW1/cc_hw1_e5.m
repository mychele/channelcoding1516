clear
disp('This is the given pcc matrix H = [B, C]perm')
H = [ 1 0 1 1 1 0 0 0 1 1;
	  0 1 0 1 1 1 1 0 0 1;
	  0 1 0 0 1 1 1 1 0 1;
	  1 1 1 1 0 0 0 1 1 0;
	  1 0 1 0 0 1 1 1 1 0; ];
disp(H);
% notice that it is full rank;
disp(strcat('Rank of H is ', num2str(rank(H))));

% build the permutation matrix
perm_matrix = eye(max(size(H)));
% in this case
row4 = perm_matrix(4, :);
row5 = perm_matrix(5, :);
row6 = perm_matrix(6, :);
perm_matrix(4, :) = row5;
perm_matrix(5, :) = row6;
perm_matrix(6, :) = row4;
disp('Permutation matrix');
disp(perm_matrix);

disp('This is Htilde=[B, C]=H*transpose(perm) with C invertible')
Htilde_nonstd = H*perm_matrix';
disp(Htilde_nonstd);

disp('Perform some row operations on Htilde to get C^{-1}')
Hb = sum_rows([eye(5), Htilde_nonstd], 1, 2)
Hb = sum_rows(Hb, 2, 3)

Hb = sum_rows(Hb, 1, 4)

Hb = sum_rows(Hb, 2, 5)
Hb = sum_rows(Hb, 4, 5)
Hb = sum_rows(Hb, 3, 4)
%Hb = sum_rows(Hb, 3, 5)


Hb = sum_rows(Hb, 5, 3)
Hb = sum_rows(Hb, 4, 3)
Hb = sum_rows(Hb, 5, 1)
Hb = sum_rows(Hb, 4, 1)
Hb = sum_rows(Hb, 4, 2)

A = Hb(:, 6:10)
C_1 = Hb(:, 1:5)
Htilde = Hb(:, 6:end)

Gtilde = [eye(5); A];
disp('This is the Gp generating matrix')
disp(Gtilde);
mod(Htilde*Gtilde, 2)
mod(Htilde_nonstd*Gtilde, 2)

disp('This is the G generating matrix')
G = perm_matrix'*Gtilde;
disp(G);
disp('and this is the result of H*G')
mod(H*G, 2)