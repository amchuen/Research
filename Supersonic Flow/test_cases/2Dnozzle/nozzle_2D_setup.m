clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
noz = load('mcCabe_nozzle.mat');
% Define Nozzle Region
options = optimset('TolX', 1e-10);
func = @(xval) ppval(noz.curve, xval);
xThroat = fminbnd(func, 0, 1, options);
x_vals = linspace(xThroat, 1, 201);
dx = x_vals(2) - x_vals(1);

% Extend x for exit region
x_vals = [x_vals, (x_vals(end)+(dx:dx:0.25))];

dy = 0.005;

y_max = max(func(x_vals));
y_vals = 0:dy:y_max;

Y_L = -func(x_vals(x_vals<=1)) ;
dyBdx = zeros(size(x_vals));
for i = 2:(length(Y_L)-1)
   dyBdx(i) = (Y_L(i+1) - Y_L(i-1))/(2*dx);
%    dyUdx(i) = (Y_U(i+1) - Y_U(i-1))/(2*dx);
end
dyBdx(1) = (-1.5*Y_L(1) + 2*Y_L(2) - 0.5*Y_L(3))./dx;
dyBdx(length(Y_L)) = (1.5*Y_L(end) - 2*Y_L(end-1) + 0.5*Y_L(end-2))./dx;



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
BC.N.range = [1, find(x_vals==x_vals(end))];

% Inlet
BC.W.varType = BC.N.varType;
BC.W.varName = BC.N.varName;
BC.W.physical = 'inlet';
BC.W.val = {1, 1, 0, 1./(FL.gam*(FL.gam-1)*FL.M0^2)+0.5};
BC.W.range = [1, find(y_vals == y_vals(end))];

% Nozzle
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.varName = BC.N.varName;
BC.S.dydx = dyBdx;
BC.S.range = [1, find(x_vals==1,1,'last')];

% % Exit Region
BC.S(2).physical = 'smallDisturb';
BC.S(2).varType = BC.N.varType;
BC.S(2).varName = BC.N.varName;
BC.S(2).perturbPrim = cat(3, zeros(size(x_vals)), -0.1*ones(size(x_vals)), repmat(zeros(size(x_vals)),1,1,2));
BC.S(2).range = [find(x_vals == 1,1,'last')+1, length(x_vals)];

% Outlet
M_e = 1.25;
p_i = p0.*((1+0.5*(FL.gam-1).*M_e^2)/(1+0.5*(FL.gam-1))).^(-FL.gam/(FL.gam-1));
BC.E.physical = 'outlet';
% BC.E.exitCond = {'p', p_i};
BC.E.varType = BC.N.varType;
BC.E.varName = BC.N.varName;
BC.E.range = [1, find(y_vals == y_vals(end))];

%% Run Simulation
GR.ratio = 1;
U0 = cat(3,rho0.*ones(size(GR.XX)), rho0.*u0.*ones(size(GR.XX)), v0.*ones(size(GR.XX)), ((rho0.^FL.gam)./(FL.gam*(FL.gam-1))+0.5.*(rho0.*u0^2)).*ones(size(GR.XX)));%cat(3, ones(size(GR.XX)), GR.XX);

% OUT = dufortFrankel(GR, FL, BC, @eulerIsenFunc, @vonNeumRichtVisc, U0);
fluxFunc = @(GR, FL, BC, EE) fluxCD_2Diff(@fullEuler, @vonNeumRichtVisc, GR, FL, BC, EE);
OUT = threeLevelExplicit(GR, FL, BC, U0, fluxFunc);

%% Post Process
close all;

BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho E'};


geomName = 'mcCabe_oblique';
folderName = ['M_' num2str(FL.M0)];
dirName = [pwd '\' geomName '\' folderName '\'];
% if ~exist(dirName, 'dir')
%     mkdir(dirName);
% end

postProcess(GR, BC, FL, OUT, dirName);
save([dirName 'OUT']);