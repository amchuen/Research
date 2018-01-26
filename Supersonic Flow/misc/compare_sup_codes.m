clc;
clear;
close all;

%% Set UP
M_inf = 1.3;

% Generate grid
grid.dx = 0.02;
grid.dy = 1.1*tan(asin(1/M_inf))*grid.dx;
grid.xvals = (-4*grid.dx:grid.dx:2);
grid.yvals = (floor(-tan(asin(1/M_inf(1)))*max(grid.xvals))-3.5*grid.dy):grid.dy:(ceil(tan(asin(1/M_inf(1)))*max(grid.xvals))+3.5*grid.dy);
% grid.yvals = 0.5*grid.dy:grid.dy:(ceil(tan(asin(1/M_inf(1)))*max(grid.xvals))+3.5*grid.dy);
grid.ii_bot = sum(grid.yvals <0);
grid.ii_top = grid.ii_bot + 1;
[grid.XX, grid.YY] = meshgrid(grid.xvals, grid.yvals);

% Create Airfoil geometry
tau = 0.1;
B_range = (grid.xvals>=0);%&(grid.xvals <=1);
xB = grid.xvals(B_range);
% yB = 2*tau.*grid.xvals((grid.xvals>=0)&(grid.xvals <=1)).*(1- grid.xvals((grid.xvals>=0)&(grid.xvals<=1))); % biconvex
yB = [0.25.*sqrt(1 - (grid.xvals((grid.xvals>=0)&(grid.xvals <= 1.0))-1).^2),...
        0.25.*ones(size(grid.xvals(grid.xvals > 1.0)))];

% Initialize Flow Parameters
flow.U_inf = 1.0;
flow.Ma = M_inf;
flow.AOA = 0.0;
arc = -1;

%% Solve 

[~, CP] = flow_supersonic(flow, grid, yB, arc);


figure();
contourf(grid.XX,grid.YY,CP, 50);
colorbar('eastoutside');