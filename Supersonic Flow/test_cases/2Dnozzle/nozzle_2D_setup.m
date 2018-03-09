clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dy = 0.0025;
dx = 0.001;

y_max = 1;
x_max = 1;%7+20*dx;
x_min = 0.5;%-7-39*dx; %(-19*dx);
x_vals = x_min:dx:x_max;
y_vals = -y_max:dy:y_max;
% x_vals = linspace(x_min, x_max, 151); 
% y_vals = linspace(-y_max, y_max, 151);

dx = x_vals(2) - x_vals(1);
dy = y_vals(2) - y_vals(1);

% Field Axis Values
[GR.XX, GR.YY] = meshgrid(x_vals, y_vals);
GR.dx = dx;
GR.dy = dy;
GR.isPolar = 0;

%% FL - fluid parameters
FL.gam = 2; % heat 
FL.M0 = 1;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-5;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = dx^2;
GR.CFL = 1;

%% Diffusion Coefficients

% epsFunc = @(FF, GR, BC, DIR) varVisc;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.1;
m_x = tand(8); % dy/dx
% x_vals = x_vals;
% dx = dx;
Y_U = 0.5.*(1 + (2.*x_vals-1).^2) + 0.5;
Y_L = -0.5.*(1 + (2.*x_vals-1).^2) + 0.5;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));
dyUdx = dyBdx;

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (Y_L(i+1) - Y_L(i-1))/(2*dx);
   dyUdx(i) = (Y_U(i+1) - Y_U(i-1))/(2*dx);
end

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'wall';
% BC.N.val = {1, 1, 0};
BC.N.varType = {'s','v2', 'v1'};
BC.N.varName = {'\rho', '\rho u', '\rho v'};
BC.N.dydx = dyUdx;

% Inlet
BC.W.varType = BC.N.varType;
BC.W.varName = BC.N.varName;
BC.W.physical = 'inlet';
BC.W.val = {1, 1, 0};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.varName = BC.N.varName;
BC.S.dydx = dyBdx;

% Outlet
BC.E.physical = 'outlet';
BC.E.exitCond = {'p', 1};
BC.E.varType = BC.N.varType;
BC.E.varName = BC.N.varName;

%% Run Simulation

U0 = cat(3,repmat(ones(size(GR.XX)),1,1,2), zeros(size(GR.XX)));%cat(3, ones(size(GR.XX)), GR.XX);

% OUT = dufortFrankel(GR, FL, BC, @eulerIsenFunc, @vonNeumRichtVisc, U0);
OUT = threeLevelExplicit(GR, FL, BC, U0, @eulerIsenFunc, @VRdiffusion);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v'};


geomName = 'nozzle2D';
folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' folderName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcess(GR, BC, FL, OUT, dirName);
save([dirName 'OUT']);