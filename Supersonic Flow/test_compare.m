clc;
close all;
clear;

%% Run tsf_airfoil.m

transonic_flow_airfoil;

close all;
clearvars -except PHI_new

%% Run tsf_cart

% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
    % Types of inputs needed
    % XX, YY
    % dx, dy
    
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% Field Axis Values
y_max = 7;
x_max = 3 + 20*GR.dx;
x_min = -39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

% % Body Values - Parabolic
tau = 0.05;
YY_B = [zeros(size(GR.x_vals(GR.x_vals <0))), ...
        2*tau.*GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1)).*(1- GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1))),...
        zeros(size(GR.x_vals(GR.x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B))
   dyBdx(i) = (YY_B(i) - YY_B(i-1))/(GR.dx);
end
    
% FL - fluid parameters
    % M0 - freestream mach number
    % gam - gamma of the flow
    
FL.gam = 1.4; % heat 
FL.M0 = 1.401;

% CT - simulation control, including tolerances, viscous factor gain, etc.
    % Types of inputs needed
    % tol - tolerances
    % v_coeff - gain for artificial viscosity
    % t_rho - add time correction to density
    % alpha - damping for three-level iteration scheme
    % CT.dt - spacing for three-level scheme
    
% CT.tol = 1e-5;
% CT.v_coeff = 1.0;
% CT.t_rho = 1.0;

CT.v_coeff = 1;
CT.tol = 1e-5;
CT.alpha = 50;
CT.dt = 0.005;

% Boundary Conditions
BC.Vy_II = zeros(size(GR.YY(end,:)));
BC.Vx_II = ones(size(GR.XX(end,:)));
BC.Vx_I = ones(size(GR.XX(:,1)));
BC.PHI_II = GR.XX(end,:);
BC.PHI_I = GR.XX(:,1);
BC.dyBdx = dyBdx;

[OUT] = tsf_cart(GR, FL, BC, CT);

