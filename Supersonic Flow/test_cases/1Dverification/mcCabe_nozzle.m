clc;
clear;
close all;

dirName = 'mcCabe_nozzle';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');
% addpath('../../Matrix Solvers/');

%% Generate Grid
noz = load('mcCabe_nozzle.mat');
xx = linspace(0,1,201);%linspace(0.5*(1-sqrt(3)/3),1,201);
dx = xx(2) - xx(1);

visc_x = 0.005;
% beta = 5e-6;
ratio_list = 0.9875; %linspace(0.975, 1, 7);% [0.85, 0.9, 0.999, 1];
tol = 1e-3;

%% Nozzle Geometry
aa = 0.6;
options = optimset('TolX', 1e-10);
func = @(xval) ppval(noz.curve, xval);
xThroat = fminbnd(func, 0, 1, options);
xx = linspace(xThroat, 1, 201);

g_x = 2.*ppval(noz.curve, xx);
dgdx = [(-1.5.*g_x(1) + 2.*g_x(2) - 0.5*g_x(3))./dx, (g_x(3:end) - g_x(1:end-2))./(2*dx), (0.5*g_x(end-2)-2.*g_x(end-1)+1.5*g_x(end))./dx];

figure();plot(xx, g_x);

%% Fluid Properties
% gam_list = 1.4:.1:2;
gam = 1.4;
cfl = 1;

%% Perform Runs
% Assign Gamma
%     gam = gam_list(ii);
ratio = 0.9875;
dt = dx*0.5;

% Initialize Field Values
% Inlet Conditions -> s0 = 0 (isentropic relations)
if xx(1) == 0
    rho0 = (0.5.*(gam+1))^(1/(gam-1));
    u0 = 0;
else %if xx(1) == 0.5
    rho0 = 1;
    u0 = 1;
end
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;
H0 = 3;

% Exit Conditions
M_e = 3;
p_i = p0.*((1+0.5*(gam-1).*M_e^2)/(1+0.5*(gam-1))).^(-gam/(gam-1));

UU = [rho0; rho0*u0; rho0*E0].*g_x;
UU = repmat(UU,1,1,3);

% Setup Time Simulation?
time = 0;
tEnd = 10*5;
[flux, Umax, ~, visc_beta] = fx_2Diff(@fluxFunc, @VRvisc, UU(:,:,end), g_x, gam, dx);
res = reshape(max(abs(flux), [], 2), 1, size(flux,1));
UU_x = UU(:,xx==0.65,3);
dtLast = dt;
beta = visc_beta.*(dt^2)/(ratio*dx^2);

figure(1);
resRho = semilogy(res(:,1)); hold on;
resU = semilogy(res(:,2));
resE = semilogy(res(:,3));
legend('\rho', 'u', 'e', 'Location', 'BestOutside');
title('Residual');
movegui(gcf, 'west');

% figure(2);
% [~, PP] = fluxFunc(UU(:,:,end)./g_x, gam);
% r_x = plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
% u_x = plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
% e_x = plot(xx, UU(3,:,3)./UU(1,:,3), 'o-');
% p_x = plot(xx, PP, '^-');
% title(['Ratio = ' num2str(ratio)]);
% legend('\rho', 'u', '\epsilon', 'P', 'Location', 'BestOutside');

while length(time) < 3 || norm(res(end,:)) > tol

   % Update Field Values and boundary conditions
    UU(:,:,1:2) = UU(:,:,2:3);

    % Check CFL
    if Umax*dt/dx ~= cfl
        dt = cfl.*dx./Umax;
        if abs(log10(dt/dtLast)) > 0.301 %1e-4
            dtLast = dt;
            fprintf('Time-step changing!\nNew time step: %0.5e\n', dt);
        end
    end

    % Calculate Next Time-Step
    beta = visc_beta.*(dt^2)/(ratio*dx^2);
    UU(:,2:end-1,3) = ((2.*beta./(dt^2)).*UU(:,2:end-1,2)-(beta./(dt^2)-0.5./dt).*UU(:,2:end-1,1)-flux)./(beta./(dt^2)+0.5/dt);
    time(end+1) = time(end)+dt;

    % Update Outflow Boundary Condition
    % 1) Extrapolate rho and E
    UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);

    % 2) Fix E
    UU(3,end,2:3) = (p_i.*g_x(end)/(gam-1) + 0.5.*(UU(2,end,2:3).^2)./UU(1,end,2:3));

    % Update Flux and Pressure
    [flux, Umax, ~, visc_beta] = fx_2Diff(@fluxFunc, @VRvisc, UU(:,:,end), g_x, gam, dx);
    res(end+1,:) = reshape(max(abs(flux), [], 2), 1, size(flux,1));

    % Plot Residuals
    resRho.YData(end+1) = res(end,1);
    resU.YData(end+1) = res(end,2);
    resE.YData(end+1) = res(end,3);
%     r_x.YData = UU(1,:,3)./g_x;
%     u_x.YData = UU(2,:,3)./UU(1,:,3);
%     e_x.YData = UU(3,:,3)./UU(1,:,3);
%     [~, PP] = fluxFunc(UU(:,:,end)./g_x, gam);
%     p_x.YData = PP;
    drawnow;

end

%% Post Process
ii=2;
fileName = ['ratio_' num2str(ii)];
save([dirName '\' fileName]);

figure(1);
saveas(gcf, [dirName '\resFig_' fileName]);

figure(2);
saveas(gcf, [dirName '\oscFig_' fileName]);

figure();
[~, PP] = fluxFunc(UU(:,:,end)./g_x, gam);
plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
plot(xx, UU(3,:,3)./UU(1,:,3), '^-');
plot(xx, PP, '.-');
title(['Ratio = ' num2str(ratio)]);
legend('\rho', 'u', 'E', 'P', 'Location', 'Best');
saveas(gcf, [dirName '\' fileName]);

fprintf('Simulation complete!\nNo. of iterations: %i\n', length(res));
if norm(res(end)) < tol
    fprintf('Success: Yes\n\n');
else
    fprintf('Success: No\n\n');
end