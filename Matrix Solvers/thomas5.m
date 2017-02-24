function xvals = thomas5(aa, bb, cc, dd, ee, ff)

% Setup
n = length(cc);
xvals = zeros(size(cc));

% Zip down
% bar denotes that first and second elements in row have been eliminated
cbar = cc(1);
dbar = dd(1);
ebar = ee(1);
fbar = ff(1);

% prime denotes that only first element in row has been eliminated
b_prime = bb(2);
c_prime = cc(2);
d_prime = dd(2);
e_prime = ee(2);
f_prime = ff(2);

for i = 2:n
    % update bars
    cbar(i) = cbar(i-1)*c_prime - b_prime * dbar(i-1);
    dbar(i) = cbar(i-1)*d_prime - b_prime * ebar(i-1);
    ebar(i) = cbar(i-1)*e_prime;
    fbar(i) = cbar(i-1)*f_prime - b_prime*fbar(i-1);
    
    if i < n
        % update primes for i+1 run
        b_prime = cbar(i-1)*bb(i+1) - aa(i+1)*dbar(i-1);
        c_prime = cbar(i-1)*cc(i+1) - aa(i+1)*ebar(i-1);
        d_prime = cbar(i-1)*dd(i+1);
        e_prime = cbar(i-1)*ee(i+1);
        f_prime = cbar(i-1)*ff(i+1) - aa(i+1)*fbar(i-1);
    end
end

% Zip up
xvals(end) = fbar(end)/cbar(end);
xvals(end-1) = (fbar(end-1) - dbar(end-1)*xvals(end))/cbar(end-1);

for ii = (n-2):-1:1
   xvals(ii) = (fbar(ii) - dbar(ii)*xvals(ii+1) - ebar(ii)*xvals(ii+2))/cbar(ii); 
end

end