clc;
clear;
close all;

dirName = 'twoThroat';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');

%% Define Grid

aftThrArea = @(x) pi.*(-sin(10.*pi.*x)./250 + x./100 + 1/100).^2;
x_vals = linspace(0.05, 0.75, 201);
g_x = aftThrArea(x_vals)./aftThrArea(0.05);
x_vals = x_vals ./aftThrArea(0.05);
dx = x_vals(2) - x_vals(1);

dgdx = [(-1.5.*g_x(1) + 2.*g_x(2) - 0.5*g_x(3))./dx, (g_x(3:end) - g_x(1:end-2))./(2*dx), (0.5*g_x(end-2)-2.*g_x(end-1)+1.5*g_x(end))./dx];

figure();plot(x_vals, g_x);
tol = 1e-6;

%% Fluid Properties
gam = 1.4;
cfl = 1;

%% Define Boundary Conditions

ratio = 0.9875;
dt = dx^2*0.5;

% Throat Conditions
rho0 = 1;
u0=1;
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5*u0^2;
H0 = (gam/(gam-1))*p0/rho0 + 0.5;

% Exit Conditions.... need to nondimensionalize the exit pressure
p_e = 0.4101;

UU = [rho0; rho0*u0; rho0*E0].*g_x;
UU = repmat(UU,1,1,3);

%% Setup Time Simulation

time = 0;
tEnd = 10*5;
[flux, Umax, ~, visc_beta] = fx_2Diff(@fluxFunc, @VRvisc, UU(:,:,end), g_x, gam, dx);
res = reshape(max(abs(flux), [], 2), 1, size(flux,1));
% UU_x = UU(:,x_vals==0.65,3);
dtLast = dt;
beta = visc_beta.*(dt^2)/(ratio*dx^2);

figure(1);
resRho = semilogy(res(:,1)); hold on;
resU = semilogy(res(:,2));
resE = semilogy(res(:,3));
legend('\rho', 'u', 'e', 'Location', 'BestOutside');
title('Residual');
movegui(gcf, 'west');

%% Run Simulation
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
%     UU(:,end,2:3) = (5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3)); % - 1/3.*UU(:,end-2,2:3);
    UU(:,end,2:3) = (2.*UU(:,end-1,2:3)./g_x(end-1) - UU(:,end-2,2:3)./g_x(end-2)).*g_x(end);

    % 2) Fix End Condition
%     p_i = (UU(1,end,2:3).*H0./(g_x(end)*gam))./(0.5*M_e^2+1/(gam-1));
%     UU(1,end,2:3) = rho_e*g_x(end);
%     UU(2,end,2:3) = rho_e*u_e*g_x(end);
%     UU(3,end,2:3) = (p_e/(gam-1) + 0.5.*rho_e.*u_e^2).*g_x(end);
%     UU(3,end,2:3) = (p_e*g_x(end)/(gam-1) + 0.5.*(UU(2,end,2:3).^2)./UU(1,end,2:3));
%     UU(2,end,2:3) = UU(1,end,2:3).*u_e;

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
% ii=3;
fileName = ['ratio_' num2str(ii)];
save([dirName '\' fileName]);

figure(1);
saveas(gcf, [dirName '\resFig_' fileName]);

% figure(2);
% saveas(gcf, [dirName '\oscFig_' fileName]);

figure();
[~, PP] = fluxFunc(UU(:,:,end)./g_x, gam);
plot(x_vals, UU(1,:,3)./g_x, '*-'); hold on;
plot(x_vals, UU(2,:,3)./UU(1,:,3), 'o-');
plot(x_vals, UU(3,:,3)./UU(1,:,3), '^-');
plot(x_vals, PP, '.-');
plot(x_vals, (UU(2,:,3)./UU(1,:,3))./sqrt(gam.*PP./(UU(1,:,3)./g_x)), '--');
title(['Exit Mach Number:' num2str(M_e)]);
legend('\rho', 'u', 'E', 'P','Ma', 'Location', 'bestoutside');
saveas(gcf, [dirName '\' fileName]);

fprintf('Simulation complete!\nNo. of iterations: %i\n', length(res));
if norm(res(end)) < tol
    fprintf('Success: Yes\n\n');
else
    fprintf('Success: No\n\n');
end