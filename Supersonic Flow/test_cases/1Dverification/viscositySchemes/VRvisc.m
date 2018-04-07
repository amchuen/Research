function epsilon = VRvisc(U, dx)

visc = 0.5;
visc_ss = 0.2*dx;
% visc_ss = 0.001;

epsilon = visc.*dx.*abs(diff(U,1,2)).^2;
denom = (0.5.*(abs(U(:,2:end,:))+abs(U(:,1:end-1,:))));
epsilon(denom~=0) = epsilon(denom~=0) ./ denom(denom~=0);

epsilon = epsilon + visc_ss;
% epsilon = visc_ss;

end