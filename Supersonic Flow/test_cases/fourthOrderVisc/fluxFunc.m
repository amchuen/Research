function flux = fluxFunc(U, dx)

% Calculate Second Derivatives
d2U = @(U) (U(3:end) - 2.*U(2:end-1) + U(1:end-2));%./(dx^2);

% Get Eigenvalue of system
C_speed = U(2:end-1);

% Calculate 2nd-Order Visc
visc2_E = VRvisc(U(2:end),dx);
visc2_W = VRvisc(U(1:end-1),dx);

% Run Checks on 2nd-Order
% if max(visc2_E(:)) > 0.5*max(abs(C_speed(:)))*dx
%     visc2_E(:) = visc2_E(:).*0.45*max(abs(C_speed(:))).*dx./(max(visc2_E(:)));
% end
% 
% if max(visc2_W(:)) > 0.5*max(abs(C_speed(:)))*dx
%     visc2_W(:) = visc2_W(:).*0.45*max(abs(C_speed(:))).*dx./(max(visc2_W(:)));
% end

% Calculate Fourth-Order Visc
incFact = 0.5;
visc4_E = (incFact.*abs(C_speed).*(dx^3)./8 - 0.25.*visc2_E.*dx^2)./(dx^3);
visc4_W = (incFact.*abs(C_speed).*(dx^3)./8 - 0.25.*visc2_W.*dx^2)./(dx^3);

visc4_E = max(0,visc4_E);
visc4_W = max(0,visc4_W);


flux =  C_speed .* (U(3:end) - U(1:end-2))./(2*dx)+...
        (visc4_E.*([d2U(U(2:end)),0] - d2U(U)) - visc4_W.*(d2U(U) - [0,d2U(U(1:end-1))]))./(dx) - ...
        (visc2_E.*(U(3:end) - U(2:end-1)) - visc2_W.*(U(2:end-1)-U(1:end-2)))./(dx^2);


end