function varargout = jst_1D(eqnFunc, UU, g_x, gam, dx) % Jameson-Schmidt-Turkel Scheme - 1 Dimensional

%% Nonlinear Difference Terms
% Get nonlinear terms from desired conservation equations
[FF, PP] = eqnFunc(UU./g_x, gam);
FF = FF.*g_x;

% Acoustic Speed (Speed of Sound)
CC = sqrt(gam.*PP.*g_x./UU(1,:));

% Calculate central difference
dgdx = (g_x(3:end) - g_x(1:end-2))./(2*dx);
F_U = (FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1;0].*PP(2:end-1).*dgdx;
    
% Check CFL
U_0 = abs(UU(2,:)./UU(1,:));
U_pA = abs(U_0 + CC);%.*g_x(2:end-1));
U_mA = abs(U_0 - CC);%.*g_x(2:end-1));
Umax = max([max(U_0(:)), max(U_pA(:)), max(U_mA(:))]);
% Umax_loc = max(max(U_0,U_pA),U_mA);

%% Calculate Viscosity 
k2 = 1/4; k4 = 1/256; c4 = 1;
Plim = 1e-7;

r_j = U_0 + CC;
r_jf = 0.5.*(r_j(3:end) + r_j(2:end-1));
r_jb = 0.5.*(r_j(2:end-1) + r_j(1:end-2));

s_j = [0, abs((PP(3:end)-2.*PP(2:end-1)+PP(1:end-2)))./max(abs(PP(3:end) - PP(2:end-1)) + abs(PP(2:end-1) - PP(1:end-2)),Plim), 0];
% s_j = [0, abs((PP(3:end)-2.*PP(2:end-1)+PP(1:end-2)))./(abs(PP(3:end)) + 2.*abs(PP(2:end-1)) + abs(PP(1:end-2))), 0];

eps2_j = k2.*r_j.*s_j;
eps2_jf = 0.5.*(eps2_j(3:end) + eps2_j(2:end-1));
eps2_jb = 0.5.*(eps2_j(1:end-2) + eps2_j(2:end-1));

eps4_jf = max(0, k4-c4.*eps2_jf./r_jf);
eps4_jb = max(0, k4-c4.*eps2_jb./r_jb);

% eps4_jf = max(0, k4.*-c4.*eps2_jf);
% eps4_jb = max(0, k4.*-c4.*eps2_jb);

%% Diffusion Terms
% Second-Order Diffusion
d2_jf = UU(:,3:end) - UU(:,2:end-1);
d2_jb = UU(:,2:end-1) - UU(:,1:end-2);

% Calculate Second-Derivatives -> assume 2nd derivs are zero at boundaries
d2Udx = [zeros(size(UU(:,1))), d2_jf - d2_jb, zeros(size(UU(:,end)))];

% fourth-order diffusion
d4_jf = d2Udx(:,3:end) - d2Udx(:,2:end-1);
d4_jb = d2Udx(:,2:end-1) - d2Udx(:,1:end-2);

%% Total Dissipation
% dt_f = dx./(0.5.*(Umax_loc(3:end) + Umax_loc(2:end-1)));
% dt_b = dx./(0.5.*(Umax_loc(2:end-1) + Umax_loc(1:end-2)));
% D_UF = (eps2_jf.*d2_jf - eps4_jf.*d4_jf); %(dx./dt_f).*
% D_UB = (eps2_jb.*d2_jb - eps4_jb.*d4_jb); %(dx./dt_b).*

D_U2 = (eps2_jf.*d2_jf - eps2_jb.*d2_jb)./dx;
D_U4 = (eps4_jf.*d4_jf - eps4_jb.*d4_jb)./dx;

D_U = D_U2 - D_U4;%(D_UF - D_UB)./dx;

%% Outputs
% Output Flux Values
flux = F_U - D_U;

varargout{1} = flux;
varargout{2} = Umax;
varargout{3} = CC;
varargout{4} = 0.5.*(eps2_jf + eps2_jb);

end