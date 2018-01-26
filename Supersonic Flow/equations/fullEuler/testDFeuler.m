clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dx = 0.04;
dy = 0.04;

% Field Axis Values
y_max = 50;
x_max = 11;%7+20*dx;
x_min = -10;%-7-39*dx; %(-19*dx);
x_vals = x_min:dx:x_max;
y_vals = 0:dy:y_max;
[GR.XX, GR.YY] = meshgrid(x_vals, y_vals);
GR.dx = dx;
GR.dy = dy;
GR.isPolar = 0;

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 0.85;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-5;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = dx^2;
GR.CFL = 1;

%% Diffusion Coefficients

epsFunc = @(GR, BC, DIR) 0.005;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.1;
m_x = tand(8); % dy/dx
% x_vals = x_vals;
% dx = dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'inlet';
BC.N.val = {1, 1, 0, 1./(FL.gam*(FL.gam-1)*FL.M0^2)+0.5};
BC.N.varType = {'s','v2', 'v1', 's'};
BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho\epsilon'};
BC.N.dydx = 0;

% Inlet
BC.W = BC.N;
% BC.W.val = {1, x_vals(1)-GR.dx};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.dydx = dyBdx;

% Outlet
BC.E.physical = 'outlet';
BC.E.varType = BC.N.varType;
% BC.E = BC.N;
% BC.E.val = {1, x_vals(end)+GR.dx};

%% Run Simulation

U0 = cat(3,repmat(ones(size(GR.XX)),1,1,2), zeros(size(GR.XX)), ones(size(GR.XX))./(FL.gam*(FL.gam-1)*FL.M0^2)+0.5);%cat(3, ones(size(GR.XX)), GR.XX);

OUT = dufortFrankel(GR, FL, BC, @eulerFullFunc, epsFunc, U0);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho e'};


geomName = 'biconvex';
folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' folderName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcess(GR, BC, FL, OUT, dirName);
save([dirName 'OUT']);