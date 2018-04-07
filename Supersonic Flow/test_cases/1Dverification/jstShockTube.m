clc;
clear;
close all;

dirName = 'jstShockTube';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');
addpath('../../../Matrix Solvers/');

%% Generate Grid

xx = linspace(0,1,1001);
dx = xx(2) - xx(1);
dt = 0.1;

visc_x = 0.0005;
% beta = 5e-6;
ratio_list = linspace(0.975, 1, 7);% [0.85, 0.9, 0.999, 1];
tol = 1e-3;

%% Nozzle Geometry
g_x = ones(size(xx));
dgdx = [0, (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];

%% Fluid Properties
gam = 1.4;
cfl = 0.95;

%% Perform Runs

ratio = 0.85;
dt = dx^2;
beta = visc_x*(dt^2)/(ratio*dx^2);

% Initialize Field Values
UU = ones([3, length(xx)]);

% density
UU(1,:) = (xx<=0.5)+0.125.*(xx>0.5); 

% velocity
UU(2,:) = 0; 

%Energy
UU(3,:) = 2.5.*(xx<=0.5)+0.25.*(xx>0.5); 

UU = repmat(UU,1,1,3);

% Setup Time Simulation?
time = 0;
tEnd = 10*5;
[flux, Umax, ~, visc_beta] = jst_1D(@fluxFunc, UU(:,:,end), g_x, gam, dx);
res = reshape(max(abs(flux), [], 2), 1, size(flux,1));
dtLast = dt;

figure(1);
resRho = semilogy(res(:,1)); hold on;
resU = semilogy(res(:,2));
resE = semilogy(res(:,3));
legend('\rho', 'u', 'Location', 'BestOutside');
title('Residual');
movegui(gcf, 'west');

figure(2);
rhox = plot(xx, UU(1,:,end)./g_x); hold on;
ux = plot(xx, UU(2,:,end)./UU(1,:,end));
ex = plot(xx, UU(3,:,end)./UU(1,:,end));
legend('\rho', 'u', 'e', 'Location', 'BestOutside');
title('Profile');
movegui(gcf, 'east');
%     Ux_fig = figure();

while length(time) < 500 && norm(res(end,:)) > tol

   % Step-forward Field Values
    UU(:,:,1:2) = UU(:,:,2:3);

    % Check CFL
    if Umax*dt/dx ~= cfl %Umax*dt/dx > 1
%             error('CFL condition not met!\n');
        dt = cfl*dx./Umax;
%         visc_beta = 0.5.*(VRvisc(UU(:,2:end,end),dx) + VRvisc(UU(:,1:end-1,end),dx));
        beta = max(abs(visc_beta(:)))*(dt^2)/(ratio*dx^2);
        if abs(log10(dt/dtLast)) > 0.301
            dtLast = dt;
            fprintf('Time-step changing!\nNew time step: %0.5e\nNew beta:%0.5e\n\n', dt, beta);
        end
    end

    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (2.*UU(:,2:end-1,2)-(1-dt*0.5./beta)*UU(:,2:end-1,1)-(dt.^2).*flux./beta)./(1+0.5*dt./beta);
%     UU(:,2:end-1,3) = (UU(:,2:end-1,1)-(2.*dt).*flux);
    time(end+1) = time(end)+dt;

    % Update Outflow Boundary Condition
%     % 1) Extrapolate rho and E
%     UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);
% 
%     % 2) Fix Rho
%     UU(1,end,2:3) = g_x(end).*((gam.*p_i).^(1/gam));
    
    % Update Flux and Pressure
    [flux, Umax, ~, visc_beta] = jst_1D(@fluxFunc, UU(:,:,end), g_x, gam, dx);
    res(end+1,:) = reshape(max(abs(flux), [], 2), 1, size(flux,1));
%     UU_x(:,end+1) = UU(:,xx==0.65,3);

    % Plot Residuals
    resRho.YData(end+1) = res(end,1);
    resU.YData(end+1) = res(end,2);
    resE.YData(end+1) = res(end,3);
    rhox.YData = UU(1,:,end)./g_x;
    ux.YData = UU(2,:,end)./UU(1,:,end);
    ex.YData = UU(3,:,end)./UU(1,:,end);
    drawnow;

end
    
%% Post Process
fileName = ['ratio_' num2str(ii)];
save([dirName '\' fileName]);

figure(1);
saveas(gcf, [dirName '\resFig_' fileName]);

figure(2);
saveas(gcf, [dirName '\oscFig_' fileName]);

figure();
[~, PP] = fluxFuncIsentropic(UU(:,:,end)./g_x, gam);
plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
plot(xx, PP, '^-');
title(['Ratio = ' num2str(ratio)]);
legend('\rho', 'u', 'P', 'Location', 'Best');
drawnow;
saveas(gcf, [dirName '\' fileName]);

fprintf('Simulation complete!\nNo. of iterations: %i\n', length(res));
if norm(res(end)) < tol
    fprintf('Success: Yes\n\n');
else
    fprintf('Success: No\n\n');
end

% close all;