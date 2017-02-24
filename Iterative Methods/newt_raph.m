function [x_n, iter, res] = newt_raph (func, dfunc, x0, tol)

% Calculate the initial guess, x1
x_n = x0 - func(x0)/dfunc(x0);
res = [abs(x_n - x0)];
iter = 1;

% Loop through until difference between guesses are less than desired
% tolerance.
while res(end) > tol
   x0 = x_n;
   x_n = x0 - func(x0)/dfunc(x0);
   iter = iter + 1;
   res(iter) = abs(x_n - x0);
   
end



end