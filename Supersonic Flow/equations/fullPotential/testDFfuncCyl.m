clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
GR.isPolar = 1;

% Define Grid
% GR.dT = 3.*pi/180*0.5;
% GR.dR = 0.133*0.25*0.5;

GR.dT = pi/180;
GR.dR = 0.133*0.5;

% Field Axis Values
R_max = 50;%GR.dR.*30+0.5;
r_cyl = 0.5;
GR.r_cyl = r_cyl;
T_max = pi;
T_min = 0.5*pi; %(-19*dx);
T_vals = T_min:GR.dT:T_max;
R_vals = r_cyl:GR.dR:R_max;
[GR.TT, GR.RR] = meshgrid(T_vals, R_vals);

GR.RR_N = [0.5.*(GR.RR(2:end,:) + GR.RR(1:(end-1),:)); 0.5.*((GR.RR(end,:)+GR.dR) + GR.RR(end,:))];
GR.RR_S = 0.5.*([2.*GR.RR(1,:)-GR.dR; GR.RR(2:end,:) + GR.RR(1:(end-1),:)]);

GR.XX = GR.RR.*cos(GR.TT);
GR.YY = GR.RR.*sin(GR.TT);
% GR.dx = dx;
% GR.dy = dy;

% figure();plot(GR.XX, GR.YY, 'b-', GR.XX', GR.YY', 'b-');axis equal

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.1;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-3;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = 1e-4;
GR.CFL = 0.25;

%% Diffusion Coefficients

epsFunc = @(GR, BC, DIR) 0.01;

%% Boundary Conditions

dyBdx = zeros(size(T_vals));

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'farfield';
% BC.N.val = {1, (GR.RR(end,:)+GR.dR + (r_cyl^2)./(GR.RR(end,:)+GR.dR)).*cos(GR.TT(end,:))};
BC.N.val = {1, (1 - (r_cyl^2)./((GR.RR(end,:)).^2)).*cos(GR.TT(end,:))};
BC.N.varType = {'s','phi'};
BC.N.varName = {'\rho','\phi'};
BC.N.dydx = 0;

% Sym
BC.E.physical = 'sym';
BC.E.varType = BC.N.varType;
BC.E.val = GR.XX(:,end);
%BC.W.val = {1, x_vals(1)-GR.dx};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.dydx = dyBdx;

% Outlet
BC.W.physical = 'outlet';
BC.W.varType = BC.N.varType;
% BC.W.physical = 'sym';
% BC.W.varType = BC.N.varType;
% BC.W.val = GR.XX(:,1);
% BC.E = BC.N;
% BC.E.val = {1, x_vals(end)+GR.dx};

%% Run Simulation

U0 = cat(3, ones(size(GR.XX)), (GR.RR + (r_cyl^2)./(GR.RR)).*cos(GR.TT));
% load('cylinder\M_1.1\testGrid2\OUT.mat','OUT');
% U0 = OUT.Uvals(:,:,:,end);
% clear OUT
OUT = dufortFrankel(GR, FL, BC, @fullPotential, epsFunc, U0);

%% Post Process
close all;

folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\cylinder\' folderName '\test5\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

OUT = postProcess(GR, BC, OUT, dirName);
save([dirName 'OUT']);
