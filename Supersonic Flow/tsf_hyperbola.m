clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% Field Axis Values
y_max = 7;
x_max = 4.5 + 20*GR.dx;
x_min = -39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.6;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 1.5;
CT.tol = 1.0e-5;
CT.alpha = 35;
CT.dt = 0.05;
CT.t_rho = 1;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.15;
m_x = tand(8); % dy/dx
x_vals = GR.x_vals;
dx = GR.dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        -tau.*sqrt(1 - 0.25.*(x_vals((x_vals>=0)&(x_vals <=2))).^2) + tau,...tau.*((0.5^3).*(x_vals((x_vals>=0)&(x_vals <2))).^3),...
        tau.*ones(size(x_vals((x_vals >2))))];
%         &(x_vals <3.5)))),...
%         zeros(size(x_vals(x_vals >=3.5)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i))/(dx);
end

BC.Vy_II = zeros(size(GR.YY(end,:)));
BC.Vx_II = ones(size(GR.XX(end,:)));
BC.Vx_I = ones(size(GR.XX(:,1)));
BC.PHI_II = GR.XX(end,:);
BC.PHI_I = GR.XX(:,1);
BC.dyBdx = dyBdx;

%% RUN CALCULATION

close all;
[OUT] = tsf_cart(GR, FL, BC, CT);

%% POST PROCESS

folderName = ['M_' num2str(FL.M0)];
geomName = 'hyperbola';

post_process(OUT, GR, geomName, folderName);