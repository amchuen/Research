clc;
clear;
close all;

%% SIM CONTROL PARAMS

% Step Sizes
dr = 0.1;
dT = 0.01*pi;
dt = 0.1*dr;

alpha = 2.5;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 10;

% Field Axis Values
r_vals = r_cyl:dr:(r_max+dr);
T_vals = 0:dT:(pi);

% Time Inputs
start = 0.0;
% stop = 100 * run.times.dt;
% dt = [0.001];


%% GRID SETUP

% Build Field
run = struct('grid',[], 'flow',[], 'times',[]);
[run.grid.TT, run.grid.RR] = meshgrid(T_vals, r_vals);
run.grid.dT = dT;
run.grid.dr = dr;

run.grid.XX = run.grid.RR .* cos(run.grid.TT);
run.grid.YY = run.grid.RR .* sin(run.grid.TT);

% Initialize flow properties
run.flow.gamma = 1.4;
run.flow.PHI_INIT = run.grid.XX; 
M0 = [0, 0.4, 0.6]; % , 0.3, 0.35

%% Loop and Solve for Convergence w/ Time

resCheck = figure();
title('Residual Check');
theta = figure();
% Initialize time grids
show_field = 1;
run.times.start = start; % seconds
for nn = 1:length(M0)
    fprintf('Simulation %i Started...\n', nn);
    if nn ~= 1
       run(nn) = run(1);
    end
    
    run(nn).flow.M0 = M0(nn);
    run(nn).flow.visc = 0.0;
%     if nn > 1
%        run(nn).flow.PHI_INIT = run(nn-1).sol.PHI(:,:,end); 
%     end
    run(nn).times.dt = dt;
    run(nn).times.alpha = alpha;%10/dt(nn);
    run(nn).times.stop = 200 * dt;
    run(nn).times.tol = 1.0e-3;

    % Run Solver
%     dbstop in flowSolve_comp if ~isreal([PHI_TT, PHI_RR, RHO_TT, RHO_RR])
    [run(nn).sol.PHI, run(nn).sol.PHI_R, run(nn).sol.PHI_T, run(nn).sol.RHO,run(nn).sol.M_IJ, run(nn).sol.res, ind, x_res, y_res] = flowSolve_comp(run(nn).grid, run(nn).flow, run(nn).times);
    fprintf('Residual: %0.5f\n', run(nn).sol.res(end));
    fprintf('Simulation %i Complete...\n', nn);
    
    % Plot Residuals
    figure(resCheck);
    hold on;
    plot(1:length(run(nn).sol.res), log10(run(nn).sol.res));
    hold off;
    
    if show_field == 1
        % Plot Fields
        PHI_X = run(nn).sol.PHI_R(:,:,end).*cos(run(nn).grid.TT) - run(nn).sol.PHI_T(:,:,end).*sin(run(nn).grid.TT)./run(nn).grid.RR;
        PHI_Y = run(nn).sol.PHI_R(:,:,end).*sin(run(nn).grid.TT) + run(nn).sol.PHI_T(:,:,end).*cos(run(nn).grid.TT)./run(nn).grid.RR;
        V2_TOT = PHI_X.^2 + PHI_Y.^2;
        Cp = 1 - V2_TOT;
        figure(); % cp plots
        contourf(run(nn).grid.XX, run(nn).grid.YY, 1-(PHI_X.^2 + PHI_Y.^2), 50); %./((RR.*cos(TT)).^2)
        title(sprintf('Pressure Coefficient Contours (Sim %i)', nn));
        colorbar('eastoutside');

        figure(); % field potential
        contourf(run(nn).grid.XX, run(nn).grid.YY, run(nn).sol.PHI(:,:,end), 50);
        title(sprintf('Field Potential (Sim %i)', nn));
        colorbar('eastoutside');
        axis equal

        figure(); % x-dir velocity plots
        contourf(run(nn).grid.XX, run(nn).grid.YY, run(nn).sol.PHI_T(:,:,end), 50); %./((RR.*cos(TT)).^2)
        title(sprintf('PHI_\theta velocity (Sim %i)', nn));
        colorbar('eastoutside');
        
        figure();
        contourf(run(nn).grid.XX, run(nn).grid.YY, run(nn).sol.M_IJ(:,:,end), 50);
        title(sprintf('Mach Number (Sim %i)', nn));
        colorbar('eastoutside');
    end
    figure(theta);
    hold on;
%     h(nn) = plot(run(nn).grid.TT(1,:), run(nn).sol.PHI_T(1,:,end));
    h(nn) = plot(run(nn).grid.TT(1,:), Cp(1,:));
    s{nn} = sprintf('Mach %0.2f', M0(nn));
    xlabel('\theta');
%     ylabel('\phi_{\theta}');
    ylabel('C_p');
    title('C_p on surface of Cylinder');
    set(gca, 'Xdir', 'reverse');
    hold off;
    
    figure();
    hold on;
    for i = 1:length(x_res)
        plot(run(nn).grid.XX(y_res(i), x_res(i)), run(nn).grid.YY(y_res(i), x_res(i)), 'o');
        hold on;
    end
    axis equal;
    title('Residual Location');
end

figure(theta);
legend(h, s);
% axis equal



%% Display Results

% PHI_X = run(nn).sol.PHI_R(:,:,end).*cos(run(nn).grid.TT) - run(nn).sol.PHI_T(:,:,end).*sin(run(nn).grid.TT);
% PHI_Y = run(nn).sol.PHI_R(:,:,end).*sin(run(nn).grid.TT) + run(nn).sol.PHI_T(:,:,end).*cos(run(nn).grid.TT);
% figure(); % cp plots
% contourf(run(nn).grid.XX, run(nn).grid.YY, 1-(PHI_X.^2 + PHI_Y.^2), 50); %./((RR.*cos(TT)).^2)
% title(sprintf('Pressure Coefficient Contours (Sim %i)', nn));
% colorbar('eastoutside');
% 
% figure(); % field potential
% contourf(run(nn).grid.XX, run(nn).grid.YY, run(nn).sol.PHI(:,:,end), 50);
% title(sprintf('Field Potential (Sim %i)', nn));
% colorbar('eastoutside');
% axis equal
% 
% figure(); % x-dir velocity plots
% contourf(run(nn).grid.XX, run(nn).grid.YY, PHI_X, 50); %./((RR.*cos(TT)).^2)
% title(sprintf('U velocity (Sim %i)', nn));
% colorbar('eastoutside');

%% Check for Convergence
