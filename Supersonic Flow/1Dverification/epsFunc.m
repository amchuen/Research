function epsilon = epsFunc(U, dx)

visc = 0.005;
visc_ss = 0.001;

epsilon = visc.*dx.*abs(diff(U,1,2)).^2;
denom = (0.5.*(U(:,2:end,:)+U(:,1:end-1,:)));
epsilon(denom~=0) = epsilon(denom~=0) ./ denom(denom~=0);

epsilon = epsilon + visc_ss;

end