function xvals = pentadiagonal(E, A, D, C, F, B)
%% Pentadiagonal Solver

% This is a pentadiagonal matrix solver. The matrix has five bands given by
% EADCF, with D being the main diagonal,E and A are the lower diagonals, 
% and C and F are the upper diagonals. B is the right-hand side solution
% The result should return x in Ax = b.

% Note that all the input vectors must be the same length. I.e....
% E is defined for i = 3:N, and zero from 1-2
% A is defined for i = 2:N, and zero for 1
% D is defined for i = 1:N,
% C is defined for i = 1:N-1, and the last element is 0
% F is defined for i = 1:N-2, and the last two elements are zero

%% Start Algorithm

% Setup
Nvals = length(D);
xvals = zeros(size(D));

% Zip Down
for i = 2:Nvals-1
    D(i) = D(i) - (A(i-1)*C(i-1)/D(i-1));
    C(i) = C(i) - (A(i-1)*F(i-1)/D(i-1));
    B(i) = B(i) - (A(i-1)*B(i-1)/D(i-1));
    
    A(i) = A(i) - (E(i-1)*C(i-1)/D(i-1));
    D(i+1) = D(i+1) - (E(i-1)*F(i-1)/D(i-1));
    B(i+1) = B(i+1) - (E(i-1)*B(i-1)/D(i-1));
end

D(Nvals) = D(Nvals) - (A(Nvals-1)*C(Nvals-1)/D(Nvals-1));
xvals(Nvals) = (B(Nvals) - (A(Nvals-1)*B(Nvals-1)/D(Nvals-1)))/D(Nvals);
xvals(Nvals-1) = (B(Nvals-1) - C(Nvals-1)*xvals(Nvals))/D(Nvals-1);

% Zip Up
for i = Nvals-2:-1:1
    xvals(i) = (B(i) - F(i)*xvals(i+2) - C(i)*xvals(i+1))/D(i);
end

end