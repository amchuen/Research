function x_n = thomas3(a, b, c, d)

% Sample code for Monday discussion - 10/2/2017

% Assume a, b, c, are the sub, main, and super diagonal of the matrix A,
% respectively. The vector d is the solution of Ax = d

n = length(b); % length of output
x_n = zeros(n,1); % initialize output array

% double letter indicates "bar" system
bb = b(1);
dd = d(1);
cc = c(1);

for i = 2:n
   bb(i) = b(i) - (a(i)*cc(i-1))/(bb(i-1));
   dd(i) = d(i) - (a(i) * dd(i-1))/(bb(i-1));
   cc(i) = c(i);
end

x_n(end) = dd(end)/bb(end);

for j = (n-1):-1:1
   x_n(j) = (dd(j) - cc(j)*x_n(j+1))/(bb(j));
end


end