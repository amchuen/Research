clc;
clear;
close all;

%% FIELD SETUP

% Step Sizes
dr = 0.1;
dT = 0.01;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 3;

% Field Axis Values
r_vals = r_cyl:dr:(r_max+dr);
T_vals = 0:dT:(pi);

% Build Field
[grid.TT, grid.RR] = meshgrid(T_vals, r_vals);
grid.dT = dT;
grid.dr = dr;

grid.XX = grid.RR .* cos(grid.TT);
grid.YY = grid.RR .* sin(grid.TT);

% Initialize flow properties
flow.gamma = 1.4;
flow.M0 = 0.0;

% Initialize time grids
times.start = 0.0; % seconds
times.dt = 0.1*dr*dT;
% times.dts = [0.1, 0.05, 0.01]*dr*dT;
times.stop = 100 * times.dt;
times.alpha = 10/(times.dt);

%% Run Solver

[PHI, PHI_R, PHI_T, RHO] = flowSolve_comp(grid, flow, times);


%% Boundary Conditions

figure();
contourf(XX, YY, 1-((PHI_R(:,:,end).^2 + PHI_T(:,:,end).^2)), 50); %./((RR.*cos(TT)).^2)
colorbar('eastoutside');

figure();
for mm = 1:length(M_vals)
    plot(T_vals, PHI_T1(mm,:));
    hold on;
end
