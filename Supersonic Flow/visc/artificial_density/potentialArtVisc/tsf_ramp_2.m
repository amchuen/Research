clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% Field Axis Values
y_max = 7;
x_max = 1 + 20*GR.dx;
x_min = -39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.401;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 1.5;
CT.tol = 1e-5;
CT.alpha = 25;
CT.dt = 0.005;
CT.t_rho = 0;   

%% Boundary Conditions

% Body Values - Ramp
tau = 0.15;
m_x = tand(25); % dy/dx
x_vals = GR.x_vals;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        m_x.*x_vals((x_vals>=0)&(x_vals <= (tau/m_x))),...
        tau.*ones(size(x_vals(x_vals > (tau/m_x))))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*GR.dx);
end

% Compare with Parabolic Airfoil
YY_P = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyPdx = zeros(size(YY_P));

for i = 2:(length(YY_P)-1)
   dyPdx(i) = (YY_P(i+1) - YY_P(i-1))/(2*GR.dx);
end

% Compare with equivalent wedge
YY_W = [zeros(size(x_vals(x_vals <0))), ...
        tau.*x_vals((x_vals>=0)&(x_vals <=0.5)),...
        tau.*(1- x_vals((x_vals>0.5)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyWdx = zeros(size(YY_W));

for i = 2:(length(YY_W))
   dyWdx(i) = (YY_W(i) - YY_W(i-1))/(GR.dx);
end

figure();plot(x_vals, YY_B, x_vals, YY_P, x_vals, YY_W);
title('Body Values');
axis equal
figure();plot(x_vals, dyBdx, x_vals, dyPdx, x_vals, dyWdx);
title('Slopes');
axis equal

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
geomName = 'ramp';

post_process(OUT, GR, geomName, folderName);