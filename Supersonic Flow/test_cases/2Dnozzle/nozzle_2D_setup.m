clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dy = 0.005;
dx = 0.005;

y_max = 0.5;
x_max = 1;%7+20*dx;
x_min = 0.5;%-7-39*dx; %(-19*dx);
x_vals = x_min:dx:x_max;
y_vals = 0:dy:y_max;

Y_U = 0.5.*(1 + (2.*x_vals-1).^2) + 0.5;
Y_L = [-0.5.*(1 + (2.*x_vals(x_vals<1)-1).^2) + 0.5, -0.5.*ones(size(x_vals(x_vals>=1)))] ;
dyBdx = zeros(size(Y_L));
% dyBdx = -(2.*x_vals(x_vals < 1)-1).*2;
% dyBdx(end) = (Y_L(end) - Y_L(end-1))./(2*dx);
% dyUdx = -dyBdx;
% dyUdx(end) = (Y_U(end) - Y_U(end-1))./(2*dx);
for i = 2:(length(Y_L)-1)
   dyBdx(i) = (Y_L(i+1) - Y_L(i-1))/(2*dx);
%    dyUdx(i) = (Y_U(i+1) - Y_U(i-1))/(2*dx);
end
dyBdx(end) = -2;

% Field Axis Values
[GR.XX, GR.YY] = meshgrid(x_vals, y_vals);
GR.dx = dx;
GR.dy = dy;
GR.isPolar = 0;

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.0;

if x_min == 0
    rho0 = 1.5;%(0.5.*(gam+1))^(1/(gam-1));
    u0 = 1/3;
    dyBdx(1) = 2;
elseif x_min == 0.5
    rho0 = 1;
    u0 = 1;
    dyBdx(1) = 0;
end
v0 = 0;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-7;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = dx^2;
GR.CFL = 1;

%% Diffusion Coefficients

% epsFunc = @(FF, GR, BC, DIR) varVisc;

%% Boundary Conditions

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'sym';
BC.N.val = GR.YY(end,:);
BC.N.varType = {'s','v2', 'v1', 's'};
BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho e'};
BC.N.dydx = 0;

% Inlet
BC.W.varType = BC.N.varType;
BC.W.varName = BC.N.varName;
BC.W.physical = 'inlet';
BC.W.val = {1, 1, 0, 1./(FL.gam*(FL.gam-1)*FL.M0^2)+0.5};

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
GR.ratio = 1;
U0 = cat(3,rho0.*ones(size(GR.XX)), rho0.*u0.*ones(size(GR.XX)), v0.*ones(size(GR.XX)), ((rho0.^FL.gam)./(FL.gam*(FL.gam-1))+0.5.*(rho0.*u0^2)).*ones(size(GR.XX)));%cat(3, ones(size(GR.XX)), GR.XX);

% OUT = dufortFrankel(GR, FL, BC, @eulerIsenFunc, @vonNeumRichtVisc, U0);
fluxFunc = @(GR, FL, BC, EE) fluxCD_2Diff(@fullEuler, @vonNeumRichtVisc, GR, FL, BC, EE);
OUT = threeLevelExplicit(GR, FL, BC, U0, fluxFunc);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho E'};


geomName = 'nozzle2D';
folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' folderName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcess(GR, BC, FL, OUT, dirName);
save([dirName 'OUT']);