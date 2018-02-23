clc;
clear;
close all;

%% Generate Grid

xx = linspace(0.5,1,51);
dx = xx(2) - xx(1);
dt = 0.001;

%% Nozzle Geometry
g_x = 1 + (2.*xx-1).^2;
% Extend nozzle
% xx(end+1:end+2) = [xx(end)+dx, xx(end)+2*dx];
% g_x(end+1:end+2) = g_x(end);

dgdx = [4*(2*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];

%% Fluid Properties
gam = 1.4;
visc_x = 0.002;
visc_t = 1e-5;

%% Initialize Field Values

% Inlet Conditions -> s0 = 0 (isentropic relations)
% rho0 = 1.5;%(0.5.*(gam+1))^(1/(gam-1));
% u0 = 1/3;
rho0 = 1;
u0 = 1;
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;

% Exit Conditions
p_i = 1;
% rho_e = 1.20135; %1.2013512490310365;
u_e = 0.416198;% 0.4162064994502009;
rho_e = 1/(g_x(end)*u_e);
E_e = p_i/((gam-1)*rho_e) + 0.5.*u_e^2;

% UU(:,end) = [rho_e; rho_e*u_e; rho_e*E0].*g_x(end);
UU = [rho0.*g_x(2:end-1), rho0.*linspace(u0, 0.4, length(g_x(2:end-1))).*g_x(2:end-1), rho0.*E0.*g_x(2:end-1)]';
UU = repmat(UU,1,3);

%% Boundary Conditions
Qbc = zeros(size(UU,1),1);
Qbc([1, length(xx)-1, 2*(length(xx)-2)+1]) = [rho0; rho0*u0; rho0*E0].*g_x(1);
Qbc([length(xx)-2, 2*(length(xx)-2), 3*(length(xx)-2)]) = [rho_e; rho_e*u_e; rho_e*E_e].*g_x(end);

Qflux = zeros(size(UU,1),1);
Qflux([1, length(xx)-1, 2*(length(xx)-2)+1]) = - g_x(1).*[rho0.*u0; 0.5*(3-gam)*rho0*u0^2 + (gam-1)*rho0*E0; -0.5*(gam-1)*rho0*u0^3+gam*rho0*u0*E0];
Qflux([length(xx)-2, 2*(length(xx)-2), 3*(length(xx)-2)]) = g_x(end).*[rho_e.*u_e; 0.5*(3-gam)*rho_e*u_e^2 + (gam-1)*rho_e*E_e; -0.5*(gam-1)*rho_e*u_e^3+gam*rho_e*u_e*E_e];

%% Initialize Matrices
d2dx = visc_x.*matDiff2Op(dx, length(xx)-2, 3);
timeRes = @(U) norm(d2dx*U(:,3) - fluxMatFunc(U(:,3), gam, dx) + wallPressure(U(:,3)./repmat(g_x(2:end-1)',3,1),dgdx(2:end-1)',gam) + visc_x.*Qbc./(dx^2) - Qflux./(2*dx));

linCoeffs = (visc_t./(dt^2) + 0.5/dt).*speye(3.*(length(xx)-2)) - d2dx;
iterRes = @(U) linCoeffs*U(:,3) + fluxMatFunc(U(:,3), gam, dx) - wallPressure(U(:,3)./repmat(g_x(2:end-1)',3,1),dgdx(2:end-1)',gam) - visc_x.*Qbc./(dx^2) + Qflux./(2*dx) - 2*visc_t.*U(:,2)./dt^2 + (visc_t./(dt^2) - 0.5/dt).*U(:,1);
jacobFunc = @(U) linCoeffs + genFluxJacob(U(:,3), gam, dx) - wallPressureJacob(U(:,3)./repmat(g_x(2:end-1)',3,1), dgdx(2:end-1)', gam);
figure(1);
plot(xx(2:end-1), UU(1:length(xx)-2));
drawnow;

%% Iterate!
tRes = timeRes(UU) ;
while  tRes > 1e-3
    
    % shift time steps
    UU(:,1:2) = UU(:,2:3);
    
    % Calculate New Time Step -> Newton's Method
    iRes = norm(iterRes(UU));
    while iRes(end) > 1e-7
        jacobian = jacobFunc(UU);
        UU(:,3) = UU(:,3) + jacobian\(-iterRes(UU));
%         UU(:,3) = UU(:,3) + bicgstab(jacobian, -iterRes(UU));
        iRes(end+1) = norm(iterRes(UU));
    end
    
    tRes(end+1) = timeRes(UU);
    
    figure(1);
    plot(xx, [Qbc(1); UU(1:length(xx)-2, 3); Qbc(length(xx)-2)]./g_x');
    drawnow;
end