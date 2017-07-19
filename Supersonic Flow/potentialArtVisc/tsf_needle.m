clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% Field Axis Values
y_max = 5;
x_max = 1.5 + 20*GR.dx;
x_min = -39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.72;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 1.5;
CT.tol = 1e-5;
CT.alpha = 20;
CT.dt = 0.01;
CT.t_rho = 1;
CT.phi_xx0 = 0;

%% Boundary Conditions

% Body Values - Needle Case
tau = 0.15; % box thickness

x_vals = GR.x_vals;
dx = GR.dx;

% Needle Dimensions
x_needle = 0.75; % length of needle
t_needle = 0.5*tau;
m_x = tand(10); % needle ramp

YY_B = [zeros(size(x_vals(x_vals <0))), ...
        m_x.*x_vals((x_vals>=0)&(x_vals <= (t_needle/m_x))),... % ramp values
        t_needle.*ones(size(x_vals((x_vals>(t_needle/m_x))&(x_vals <= x_needle)))),...
        tau.*sqrt(1 - 4.*(x_vals((x_vals>x_needle)&(x_vals <= (0.5+x_needle)))-(0.5+x_needle)).^2),...
        tau.*ones(size(x_vals(x_vals > (0.5+x_needle))))];

dyBdx = zeros(size(YY_B));
for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i) - YY_B(i-1))/(dx);
end

BC.Vy_II = zeros(size(GR.YY(end,:)));
BC.Vx_II = ones(size(GR.XX(end,:)));
BC.Vx_I = ones(size(GR.XX(:,1)));
BC.PHI_II = GR.XX(end,:);
BC.PHI_I = GR.XX(:,1);
BC.dyBdx = dyBdx;

%% RUN CALCULATION

[OUT] = tsf_cart(GR, FL, BC, CT);

%% POST PROCESS

folderName = ['M_' num2str(FL.M0)];
geomName = 'needle';

post_process(OUT, GR, geomName, folderName);