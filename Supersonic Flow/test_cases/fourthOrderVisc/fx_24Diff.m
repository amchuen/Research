function varargout = fx_24Diff(eqnFunc, viscFunc, UU, g_x, gam, dx)

%% Nonlinear Difference Terms
% Get nonlinear terms from desired conservation equations
[FF, PP, CC] = eqnFunc(UU./g_x, gam);
FF = FF.*g_x;

% Calculate central difference
dgdx = [(-1.5.*g_x(1) + 2.*g_x(2) - 0.5*g_x(3))./dx, (g_x(3:end) - g_x(1:end-2))./(2*dx), (0.5*g_x(end-2)-2.*g_x(end-1)+1.5*g_x(end))./dx];
F_U = (FF(:,3:end) - FF(:,1:end-2))./(2.*dx);
F_U(2,:) = F_U(2,:) - PP(2:end-1).*dgdx(2:end-1);

% Get Eigenvalue of system
C_speed = max(CC,[],1);

%% Dissipation Terms

% Calculate Second Derivatives
d2U = [zeros(size(UU(:,1))), (UU(:,3:end) - 2.*UU(:,2:end-1) + UU(:,1:end-2)), zeros(size(UU(:,end)))];%./(dx^2);

% Calculate 2nd-Order Visc
visc2_E = viscFunc(UU(:,2:end),dx);%.*abs(max(C_speed(2:end-1), C_speed(3:end)));
visc2_W = viscFunc(UU(:,1:end-1),dx);%.*abs(max(C_speed(2:end-1), C_speed(1:end-2)));

% Run Checks on 2nd-Order
% if max(visc2_E(:)) > 0.5*max(abs(C_speed(:)))*dx
%     visc2_E(:) = visc2_E(:).*0.45*max(abs(C_speed(:))).*dx./(max(visc2_E(:)));
% end
% 
% if max(visc2_W(:)) > 0.5*max(abs(C_speed(:)))*dx
%     visc2_W(:) = visc2_W(:).*0.45*max(abs(C_speed(:))).*dx./(max(visc2_W(:)));
% end

% Calculate Fourth-Order Visc
incFact = 1/4;
visc4_E = (incFact.*abs(max(C_speed(2:end-1), C_speed(3:end)))./8 - 0.25.*visc2_E./dx);
visc4_W = (incFact.*abs(max(C_speed(2:end-1), C_speed(1:end-2)))./8 - 0.25.*visc2_W./dx);

visc4_E = max(0,visc4_E);
visc4_W = max(0,visc4_W);


flux =  F_U+...
        (visc4_E.*(d2U(:,3:end) - d2U(:,2:end-1)) - visc4_W.*(d2U(:,2:end-1) - d2U(:,1:end-2)))./(dx) - ...
        (visc2_E.*(UU(:,3:end) - UU(:,2:end-1)) - visc2_W.*(UU(:,2:end-1)-UU(:,1:end-2)))./(dx^2);

%% Output
varargout{1} = flux;
varargout{2} = max(C_speed(:));
varargout{3} = CC;
varargout{4} = 0.5.*(visc2_E + visc2_W);


end