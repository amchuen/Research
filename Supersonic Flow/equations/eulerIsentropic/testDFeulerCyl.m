clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
GR.isPolar = 1;

% Define Grid
GR.dT = 3.*pi/180*0.5;
GR.dR = 0.133*0.25*0.5;

% Field Axis Values
R_max = 25;%GR.dR.*30+0.5;
r_cyl = 0.5;
T_max = pi;
T_min = 0.5*pi; %(-19*dx);
T_vals = T_min:GR.dT:T_max;
R_vals = (r_cyl+0.5*GR.dR):GR.dR:R_max;
[GR.TT, GR.RR] = meshgrid(T_vals, R_vals);

GR.RR_N = [0.5.*(GR.RR(2:end,:) + GR.RR(1:(end-1),:)); 0.5.*((GR.RR(end,:)+GR.dR) + GR.RR(end,:))];
GR.RR_S = 0.5.*([2.*GR.RR(1,:)-GR.dR; GR.RR(2:end,:) + GR.RR(1:(end-1),:)]);
GR.r_cyl = r_cyl;

GR.XX = GR.RR.*cos(GR.TT);
GR.YY = GR.RR.*sin(GR.TT);
% GR.dx = dx;
% GR.dy = dy;

%% FL - fluid parameters
FL.gam = 2; % heat 
% FL.M0 = 1.1;
FL.M0 = 1.4;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-4;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = 0.1;
GR.CFL = 1;

%% Diffusion Coefficients

epsFunc = @(GR, BC, DIR) 0.01;

%% Boundary Conditions

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'inlet';
BC.N.val = {1,...
            (1 - (r_cyl^2)./((GR.RR(end,:)+GR.dR).^2)).*cos(GR.TT(end,:)),...
            -(1+(r_cyl^2)./((GR.RR(end,:)+GR.dR).^2)).*sin(GR.TT(end,:))};
BC.N.varType = {'s','v1', 'v2'};
BC.N.varName = {'\rho', '\rho u_r', '\rho v_\theta'};
BC.N.dydx = 0;

% Outlet
BC.W.physical = 'outlet';
% BC.W.physical = 'outlet';
BC.W.varType = BC.N.varType;
BC.W.val = GR.XX(:,1);
% BC.W.val = {1, x_vals(1)-GR.dx};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.dydx = zeros(size(GR.TT(1,:)));

% Outlet
BC.E.physical = 'sym';
BC.E.varType = BC.N.varType;
BC.E.val = GR.XX(:,end);
% BC.E = BC.N;
% BC.E.val = {1, x_vals(end)+GR.dx};

%% Run Simulation

U0 = cat(3,ones(size(GR.XX)), (1-(r_cyl./GR.RR).^2).*cos(GR.TT), -(1+(r_cyl./GR.RR).^2).*sin(GR.TT));%cat(3, ones(size(GR.XX)), GR.XX);

OUT = dufortFrankel(GR, FL, BC, @eulerIsenFunc, epsFunc, U0);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v'};


geomName = 'cylinder';
folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' folderName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcess(GR, BC, FL, OUT, dirName);
save([dirName 'OUT']);