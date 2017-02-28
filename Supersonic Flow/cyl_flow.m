clc;
clear;
close all;

%% SIM CONTROL PARAMS

% Step Sizes
dr = 0.01;
dT = 0.01;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 3;

% Field Axis Values
r_vals = r_cyl:dr:(r_max+dr);
T_vals = 0:dT:(pi);

% Time Inputs
start = 0.0;
% stop = 100 * run.times.dt;
dt = [0.1]*dr*dT;


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
run.flow.M0 = 0.0;

%% Loop and Solve for Convergenc w/ Time

resCheck = figure();
% Initialize time grids
run.times.start = start; % seconds
for nn = 1:length(dt)
    fprintf('Simulation %i Started...\n', nn);
    if nn ~= 1
       run(nn) = run(1);
    end
    
    run(nn).times.dt = dt(nn);
    run(nn).times.alpha = 10/dt(nn);
    run(nn).times.stop = 100 * dt(nn);

    % Run Solver
    [run(nn).sol.PHI, run(nn).sol.PHI_R, run(nn).sol.PHI_T, run(nn).sol.RHO, run(nn).sol.res] = flowSolve_comp(run(nn).grid, run(nn).flow, run(nn).times);
    fprintf('Residual: %0.5f\n', log(run(nn).sol.res(end)));
    fprintf('Simulation %i Complete...\n', nn);
    
    % Plot Residuals
    figure(resCheck);
    hold on;
    loglog(0:1:(100-3), run(nn).sol.res);
    hold off;
    
    % Plot Fields
    PHI_X = run(nn).sol.PHI_R(:,:,end).*cos(run(nn).grid.TT) - run(nn).sol.PHI_T(:,:,end).*sin(run(nn).grid.TT);
    PHI_Y = run(nn).sol.PHI_R(:,:,end).*sin(run(nn).grid.TT) + run(nn).sol.PHI_T(:,:,end).*cos(run(nn).grid.TT);
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
    contourf(run(nn).grid.XX, run(nn).grid.YY, PHI_X, 50); %./((RR.*cos(TT)).^2)
    title(sprintf('U velocity (Sim %i)', nn));
    colorbar('eastoutside');
    
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
