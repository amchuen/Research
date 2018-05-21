clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
noz = load('tripleThroat.mat');
% Define Nozzle Region
xThroat = 0.05;
xEnd = 0.55;
channel_width = noz.aftThrArea(xThroat);
func = @(x) 0.5.*noz.aftThrArea(x)./channel_width;
x_vals = linspace(xThroat, xEnd, 9001)./channel_width;
dx = x_vals(2) - x_vals(1);

y_max = max(func(x_vals.*channel_width));
y_vals = linspace(0,y_max, 11);
dy = y_vals(2) - y_vals(1);

g_x = -func(x_vals.*channel_width) ;
dyBdx = [(-1.5.*g_x(1) + 2.*g_x(2) - 0.5*g_x(3))./dx, (g_x(3:end) - g_x(1:end-2))./(2*dx), (0.5*g_x(end-2)-2.*g_x(end-1)+1.5*g_x(end))./dx];

% Field Axis Values
[GR.XX, GR.YY] = meshgrid(x_vals, y_vals);
GR.dx = dx;
GR.dy = dy;
GR.isPolar = 0;

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.0;

rho0 = 1;
u0 = 1;
v0 = 0;
p0 = (rho0^FL.gam)/FL.gam;
E0 = p0/((FL.gam-1)*rho0) + 0.5*u0^2;

%% Simulation control, including tolerances, viscous factor gain, etc.

GR.tol = 1e-5;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = dx^2;
GR.CFL = 1;

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
BC.N.range = [1, find(x_vals==x_vals(end))];

% Inlet
BC.W.varType = BC.N.varType;
BC.W.varName = BC.N.varName;
BC.W.physical = 'inlet';
BC.W.val = {rho0, rho0*u0, rho0*v0, rho0.*E0};
BC.W.range = [1, find(y_vals == y_vals(end))];

% Nozzle
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.varName = BC.N.varName;
BC.S.dydx = dyBdx;
BC.S.range = [1, length(x_vals)];

% % % Exit Region
% BC.S(2).physical = 'smallDisturb';
% BC.S(2).varType = BC.N.varType;
% BC.S(2).varName = BC.N.varName;
% BC.S(2).perturbPrim = cat(3, zeros(size(x_vals)), -1.5.*ones(size(x_vals)), repmat(zeros(size(x_vals)),1,1,2));
% BC.S(2).range = [find(x_vals == 1,1,'last')+1, length(x_vals)];

% AA = 0.5; CC = -3; BB = (FL.gam/(FL.gam-1))*noz.p_e_ratio.*func(xEnd)./func(xThroat);
% u_E = (-BB + sqrt(BB^2-4*AA*CC))/(2*AA);

% Outlet
BC.E.physical = 'outlet';
% BC.E.exitCond = {'p', noz.p_e_ratio};
BC.E.exitCond = {'u', noz.u_1};
BC.E.varType = BC.N.varType;
BC.E.varName = BC.N.varName;
BC.E.range = [1, find(y_vals == y_vals(end))];

%% Run Simulation
GR.ratio = 1;
% oldResults = load('nonUnique_nozzle\M_1\OUT.mat');
% U0 = oldResults.OUT.Uvals(:,:,:,end);
U0 = cat(3, rho0.*ones(size(GR.XX)), 1.75*rho0.*u0.*ones(size(GR.XX)), v0.*ones(size(GR.XX)), ((rho0.^FL.gam)./(FL.gam*(FL.gam-1))+0.5.*(rho0.*(u0)^2)).*ones(size(GR.XX)));%cat(3, ones(size(GR.XX)), GR.XX);

fluxFunc = @(GR, FL, BC, EE) fluxCD_2Diff(@fullEuler, @vonNeumRichtVisc, GR, FL, BC, EE);
OUT = threeLevelExplicit(GR, FL, BC, U0, fluxFunc);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho E'};

geomName = 'three_throat';
% geomName = 'twoShock';
resultName = 'oneExitShock';
% folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' resultName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcessHalfNozzle(GR, BC, FL, OUT, dirName, @fullEuler);
save([dirName 'OUT']);