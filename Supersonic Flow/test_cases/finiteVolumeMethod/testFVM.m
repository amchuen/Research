clc;
clear;
close all;

%% Setup

% Grid
GR = load('quad_nozzle.mat');

% Fluid Parameters
FL.M0 = 1.0;
FL.gam = 1.4;

% Control Params
GR.tol = 1e-5;
GR.tEnd = 0.2; % 10 seconds maximum?
GR.dt = 5e-3;
GR.CFL = 1;%2.*sqrt(2);

% Initialize Values
% Inlet Conditions -> s0 = 0 (isentropic relations)
rho0 = 1;
u0 = 1;
p0 = (rho0^FL.gam)/FL.gam;
E0 = p0/((FL.gam-1)*rho0) + 0.5.*u0^2;

% three-level time
UU = repmat(cat(3,rho0, rho0*u0, 0,rho0*E0).*ones(size(GR.UU_1(:,:,1))),1,1,1,3).*GR.UU_vol;

% Boundary Conditions
% Far-field... need to update this?
BC.N.physical = 'inlet';
BC.N.val = {1, 1, 0, 1./(FL.gam*(FL.gam-1)*FL.M0^2)+0.5};
BC.N.varType = {'s','v2', 'v1', 's'};
% BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho\epsilon'};
BC.N.varName = {'\rho', '\rho u', '\rho v', '\rho e'};
BC.N.dydx = 0;

% Inlet
BC.W = BC.N;
% BC.W.val = {1, x_vals(1)-GR.dx};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
% BC.S.dydx = dyBdx;

% Outlet
BC.E.physical = 'outlet';
BC.E.varType = BC.N.varType;

%% Start Simulation

res = ones(1,size(UU,3));

while norm(res(end,:)) > GR.tol
    % Update Values
    UU(:,:,:,1:2) = UU(:,:,:,2:3);

    % Update Boundaries


    % Calculate Flux Values
    [FF, GG, PP] = fullEuler(FL, BC, UU(:,:,:,2));
    
    FF_N =  GR.wt_N(:,:,1).*[FF(1:end-1,:,:); zeros(size(FF(end,:,:)))] +...
            GR.wt_N(:,:,2).*[FF(2:end,:,:); zeros(size(FF(end,:,:)))];
    GG_N =  GR.wt_N(:,:,1).*[GG(1:end-1,:,:); zeros(size(FF(end,:,:)))] +...
            GR.wt_N(:,:,2).*[GG(2:end,:,:); zeros(size(FF(end,:,:)))];
    
    FF_S =  GR.wt_S(:,:,1).*[zeros(size(FF(1,:,:))); FF(2:end,:,:)] +...
            GR.wt_S(:,:,2).*[zeros(size(FF(1,:,:))); FF(1:end-1,:,:)];
    GG_S =  GR.wt_S(:,:,1).*[zeros(size(GG(1,:,:))); GG(2:end,:,:)] +...
            GR.wt_S(:,:,2).*[zeros(size(GG(1,:,:))); GG(1:end-1,:,:)];
%     flux_N = 
    
    
end