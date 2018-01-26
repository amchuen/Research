function xvals = thomas3(aa, bb, cc, dd)
% This function takes in one array of main diagonal elements (b), and
% two sub-diagonal elements above (c) and below (a) the main diagonal.
% The final variable contains the vector B in Ax=B. Using the Thomas
% algorithm, this function returns the solutions to a tridiagonal system.

% Setup
n = length(aa);
xvals = zeros(size(aa));

% Zip up
cbar = [cc(1)/bb(1)];
dbar = [dd(1)/bb(1)];

for i = 2:(n)
    cbar(i) = cc(i)/(bb(i) - aa(i) * cbar(i-1));
    dbar(i) = (dd(i) - aa(i)*dbar(i-1))/(bb(i) - aa(i)*cbar(i-1));
end

% Zip down
xvals(end) = dbar(end);

for ii = linspace(n-1,1, n-1)
   xvals(ii) = dbar(ii) - cbar(ii)*xvals(ii+1); 
end

end