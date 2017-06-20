clc;
clear;
close all;

%% Initialize Variables

M_inf = 1.4;

% Generate grid
grid.dx = 0.02;
grid.dy = 1.1*tan(asin(1/M_inf(1)))*grid.dx;
grid.xvals = (-4*grid.dx:grid.dx:2);
grid.yvals = (floor(-tan(asin(1/M_inf(1)))*max(grid.xvals))-3.5*grid.dy):grid.dy:(ceil(tan(asin(1/M_inf(1)))*max(grid.xvals))+3.5*grid.dy);
% grid.yvals = 0.5*grid.dy:grid.dy:(ceil(tan(asin(1/M_inf(1)))*max(grid.xvals))+3.5*grid.dy);
grid.ii_bot = sum(grid.yvals <0);
grid.ii_top = grid.ii_bot + 1;
[grid.XX, grid.YY] = meshgrid(grid.xvals, grid.yvals);

% Create Airfoil geometry
tau = 0.1;
B_range = (grid.xvals>=0)&(grid.xvals <=1);
xB = grid.xvals(B_range);
yB = 2*tau.*grid.xvals((grid.xvals>=0)&(grid.xvals <=1)).*(1- grid.xvals((grid.xvals>=0)&(grid.xvals<=1))); % biconvex