clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dr = 0.05;
GR.dT = 0.01*pi;

% Field Axis Values
r_cyl = 1.0;
r_max = 8.0;

% Field Axis Values
GR.r_vals = (r_cyl+0.5*GR.dr):GR.dr:(r_max);

%     T_vals = 0:dT:(pi);  
GR.T_vals = 0:GR.dT:(pi);  

% Grid Setup
[GR.TT, GR.RR] = meshgrid(GR.T_vals, GR.r_vals);
GR.XX = GR.RR .* cos(GR.TT);
GR.YY = GR.RR .* sin(GR.TT);
% [GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 0.53;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 0.001;
CT.tol = 1.0e-5;
CT.alpha = 35;
CT.dt = 5;
CT.t_rho = 1;

%% Boundary Conditions

BC.Vr_II = cos(GR.TT(end,:)).*(1 - (r_cyl^2)./(GR.RR(end,:).^2));
BC.dirichlet.PHI_II = (GR.RR(end,:) + (r_cyl^2)./(GR.RR(end,:))).*cos(GR.TT(end,:));
visc_on = 0;

%% RUN CALCULATION

close all;
[OUT] = tsf_euler(GR, FL, BC, CT);

%% POST PROCESS

folderName = ['M_' num2str(FL.M0)];
geomName = 'euler_cyl';

post_process(OUT, GR, geomName, folderName);