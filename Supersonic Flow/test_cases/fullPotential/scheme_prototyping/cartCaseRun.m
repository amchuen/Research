clc;
close all;
clear;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-5;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = 5e-5;
GR.CFL = 1;
GR.isPolar = 0;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dx = 0.05;
dy = 0.08;

% Grid Size -> Subsonic Domain
y_max = 50;
x_max = 11;%7+20*dx;
x_min = -10;%-7-39*dx; %(-19*dx);
x_vals = x_min:dx:x_max;
y_vals = 0:dy:y_max;
[GR.XX, GR.YY] = meshgrid(x_vals, y_vals);
GR.dx = dx;
GR.dy = dy;

% % Supersonic Domain
% y_max = 30;
% x_min = -5;
% x_max = 11;
% x_vals = x_min:dx:x_max;
% y_vals = 0:dy:y_max;
% GRsup = GR;
% [GRsup.XX, GRsup.YY] = meshgrid(x_vals, y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat
M0 = [0.85, 0.9, 0.95, 1.1, 1.4, 2];
% FL.M0 = 0.98;

%% Diffusion Coefficients

epsFunc = @(GR, BC, DIR) 0.01;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.1;
% x_vals = x_vals;
% dx = dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

% Diamond Airfoil
YY_D = [zeros(size(x_vals(x_vals <0))), ...
        tau.*x_vals((x_vals>=0)&(x_vals <0.5)),...
        tau.*(1- x_vals((x_vals>=0.5)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyDdx = zeros(size(YY_D));

for i = 2:(length(YY_D)-1)
   dyDdx(i) = (YY_D(i+1) - YY_D(i-1))/(2*dx);
end

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC_B.N.physical = 'inlet';
BC_B.N.val = {1, GR.XX(end,:)};
BC_B.N.varType = {'s','phi'};
BC_B.N.varName = {'\rho','\phi'};
BC_B.N.dydx = 0;

% Inlet
BC_B.W = BC_B.N;
BC_B.W.val = {1, GR.XX(1,1)-GR.dx};

% Wall
BC_B.S.physical = 'wall';
BC_B.S.varType = BC_B.N.varType;
BC_B.S.dydx = dyBdx;

% Outlet
% BC.E.physical = 'outlet';
% BC.E.varType = BC.N.varType;
BC_B.E = BC_B.N;
BC_B.E.val = {1, GR.XX(1,end)+GR.dx};

% Diamond Airfoil Boundary Conditions
BC_D = BC_B;
BC_D.S.dydx = dyDdx;

%% Run Simulation and Post Process
 
U0 = cat(3, ones(size(GR.XX)), GR.XX);

for i = 4:length(M0)
    FL.M0 = M0(i);
    
    if M0(i) < 1
%         GR = GR;
%         BC_B.N.val{2} = GR.XX(end,:);
        BC_B.E.physical = 'inlet';
%         BC_B.E.val{2} = GR.XX(1,end)+GR.dx;
%         BC_B.W.val{2} = GR.XX(1,1)-GR.dx;
        BC_D.E.physical = 'inlet';
%         BC_D.E.val{2} = BC_B.E.val{2};
%         BC_D.W.val{2} = BC_B.W.val{2};
    else 
%         GR = GRsup;
%         BC_B.N.val{2} = GR.XX(end,:);
        BC_B.E.physical = 'outlet';
%         BC_B.W.val{2} = GR.XX(1,1)-GR.dx;
        BC_D.E.physical = 'outlet';
%         BC_D.W.val{2} = BC_B.W.val{2};
    end
    
    % Run Biconvex Case
    OUT_B = dufortFrankel(GR, FL, BC_B, @CIPM, epsFunc, U0);
    close all;
    
    % Run Diamond Case
    OUT_D = dufortFrankel(GR, FL, BC_D, @CIPM, epsFunc, U0);
    close all;

    folderName = ['M_' num2str(FL.M0)];
    dirName_B = [pwd '\biconvex\' folderName '\'];
    dirName_D = [pwd '\diamond\' folderName '\'];
    
    % Post Process

    postProcess(GR, BC_B, FL, OUT_B, dirName_B);
    close all;
    OUT = OUT_B;
    dirName = dirName_B;
    save([dirName_B 'results'], 'OUT', 'GR', 'BC_B', 'FL', 'epsFunc', 'YY_B', 'dirName');
    
    postProcess(GR, BC_D, FL, OUT_D, dirName_D);
    close all;
    OUT = OUT_D;
    dirName = dirName_D;
    save([dirName_B 'results'], 'OUT', 'GR', 'BC_D', 'FL', 'epsFunc', 'YY_D', 'dirName');
    
end