function varargout = fx_2Diff(eqnFunc, viscFunc, UU, g_x, gam, dx) % flux 2nd order diffusion

%% Nonlinear Difference Terms
% Get nonlinear terms from desired conservation equations
[FF, PP] = eqnFunc(UU./g_x, gam);
FF = FF.*g_x;

% Acoustic Speed (Speed of Sound)
CC = sqrt(gam.*PP.*g_x./UU(1,:));

% Calculate central difference
dgdx = (g_x(3:end) - g_x(1:end-2))./(2*dx);
F_U = (FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1].*PP(2:end-1).*dgdx;
    
% Check CFL
U_0 = abs(UU(2,:)./UU(1,:));
U_pA = abs(U_0 + CC);%.*g_x(2:end-1));
U_mA = abs(U_0 - CC);%.*g_x(2:end-1));
Umax = max([max(U_0(:)), max(U_pA(:)), max(U_mA(:))]);

%% 2nd Order Dissipation
% Calculate Viscosities
visc_E = viscFunc(UU(:,2:end),dx);
visc_W = viscFunc(UU(:,1:end-1),dx);

% Calculate 2nd Order Diffusion
D_U = (visc_E.*(UU(:,3:end) - UU(:,2:end-1)) - visc_W.*(UU(:,2:end-1) - UU(:,1:end-2)))./(dx.^2);

%% Outputs
% Output Flux Values
flux = F_U - D_U;

varargout{1} = flux;
varargout{2} = Umax;
varargout{3} = CC;

end