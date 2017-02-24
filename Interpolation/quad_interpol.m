function [xout, yout] = quad_interpol(xin, yin, k_init, res)

kvals = [k_init];
hvals = diff(xin); % get the various dx for the set of points

% build array of k values
for i = 1:(length(xin)-1)
    kvals(i+1) = 2 * (yin(i+1) - yin(i))/hvals(i) - kvals(i);
end

% initialize output
xout = [];
yout = [];
for ii = 1:(length(xin)-1)
   % get domain between xin(i) and xin(i+1)
   xvals = linspace(xin(ii), xin(ii+1), res+1);
   
   % solve for important constants
   aa = yin(ii);
   bb = kvals(ii);
   cc = 0.5 * (kvals(ii+1) - kvals(ii))/ hvals(ii);
   
   % get local equation and calculate values
   yvals = aa + bb.*(xvals - xin(ii)) + cc.*(xvals - xin(ii)).^2; 
   
   % build output values
   xout = [xout xvals];
   yout = [yout yvals];
end


end