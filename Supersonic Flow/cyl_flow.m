clc;
clear;
close all;

%% SIM CONTROL PARAMS

% Step Sizes
dr = 0.1;
dT = 0.01*pi;
dt = 0.1*dr;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 3;

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
M0 = [0.0, 0.2, 0.3, 0.35];

%% Loop and Solve for Convergence w/ Time

resCheck = figure();
title('Residual Check');
theta = figure();
% Initialize time grids
run.times.start = start; % seconds
for nn = 1:length(M0)
    fprintf('Simulation %i Started...\n', nn);
    if nn ~= 1
       run(nn) = run(1);
    end
    
    run(nn).flow.M0 = M0(nn);
    if nn > 1
       run(nn).flow.PHI_INIT = run(nn-1).sol.PHI(:,:,end); 
    end
    run(nn).times.dt = dt;
    run(nn).times.alpha = 5*1.5;%10/dt(nn);
    run(nn).times.stop = 202 * dt;

    % Run Solver
    [run(nn).sol.PHI, run(nn).sol.PHI_R, run(nn).sol.PHI_T, run(nn).sol.RHO, run(nn).sol.res] = flowSolve_comp(run(nn).grid, run(nn).flow, run(nn).times);
    fprintf('Residual: %0.5f\n', log(run(nn).sol.res(end)));
    fprintf('Simulation %i Complete...\n', nn);
    
    % Plot Residuals
    figure(resCheck);
    hold on;
    plot((0:1:((run(nn).times.stop - run(nn).times.start)/run(nn).times.dt - 3)), log10(run(nn).sol.res));
    hold off;
    
    % Plot Fields
    PHI_X = run(nn).sol.PHI_R(:,:,end).*cos(run(nn).grid.TT) - run(nn).sol.PHI_T(:,:,end).*sin(run(nn).grid.TT);
    PHI_Y = run(nn).sol.PHI_R(:,:,end).*sin(run(nn).grid.TT) + run(nn).sol.PHI_T(:,:,end).*cos(run(nn).grid.TT);
%     figure(); % cp plots
%     contourf(run(nn).grid.XX, run(nn).grid.YY, 1-(PHI_X.^2 + PHI_Y.^2), 50); %./((RR.*cos(TT)).^2)
%     title(sprintf('Pressure Coefficient Contours (Sim %i)', nn));
%     colorbar('eastoutside');
    
%     figure(); % field potential
%     contourf(run(nn).grid.XX, run(nn).grid.YY, run(nn).sol.PHI(:,:,end), 50);
%     title(sprintf('Field Potential (Sim %i)', nn));
%     colorbar('eastoutside');
%     axis equal
    
%     figure(); % x-dir velocity plots
%     contourf(run(nn).grid.XX, run(nn).grid.YY, PHI_X, 50); %./((RR.*cos(TT)).^2)
%     title(sprintf('U velocity (Sim %i)', nn));
%     colorbar('eastoutside');
    figure(theta);
    hold on;
    plot(run(nn).grid.TT(1,:), run(nn).sol.PHI_T(1,:,end));
    hold off;
end




%% Display Results

% figure();
% contourf(grid.XX, grid.YY, 1-((PHI_R(:,:,end-1).^2 + PHI_T(:,:,end-1).^2)), 50); %./((RR.*cos(TT)).^2)
% colorbar('eastoutside');
% 
% figure();
% contourf(grid.XX, grid.YY, PHI_R(:,:,end-1).*cos(grid.TT) - PHI_T(:,:,end-1).*sin(grid.TT), 50); %./((RR.*cos(TT)).^2)
% colorbar('eastoutside');
% 
% figure();
% quiver(grid.XX, grid.YY, PHI_R(:,:,end-1).*cos(grid.TT) - PHI_T(:,:,end-1).*sin(grid.TT), PHI_R(:,:,end-1).*sin(grid.TT) + PHI_T(:,:,end-1).*cos(grid.TT));

%% Check for Convergence
