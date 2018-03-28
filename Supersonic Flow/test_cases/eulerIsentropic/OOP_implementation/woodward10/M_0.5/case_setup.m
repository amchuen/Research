clc;
clear;
close all;

case_name = 'woodward';

%% Define grid_lims

dx = 0.05;
dy = 0.08;

% Field Axis Values - body fitted grid
grid_lims = [   0,  10, dy;
                -39*dx, 7+20*dx,    dx];... % y-range

%% Airfoil Geometry

tau = 0.1;
x_vals = grid_lims(2,1):grid_lims(2,3):grid_lims(2,2);
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i+1) = (YY_B(i) - YY_B(i-1))/(2*dx);
end

%% Control Params - simulation control, including tolerances, viscous factor gain, etc.

CT.eps_s = 0.07525; % spatial diffusion term
% eps_t = 0.005; % time diffusion term
CT.tol = 1e-5;
CT.dt = 0.002;
CT.iter_min = 300;
CT.CFL_on = 1;
CT.use_1visc = 1;
CT.is_polar = 0;
CT.case_name = 'cylinder_vec_1visc';

%% Fluid Params
FL.M0 = 0.5;
FL.gam = 1.4;

%% Boundary Condition Setup

BC_setup = {'W',        'E',        'N',            'S';...
            'inlet',    'outlet',   'far-field',    'wall';...
            1,          0,          1,              0;...
            1,          0,          1,              0;...
            0,          0,          0,              (dyBdx)};
        
test = eulerIsentropicField(CT, grid_lims, FL, BC_setup);

test = test.timeStep_tl;

%% Post Process
folderName = ['M_' num2str(FL.M0)];
geomName = [case_name num2str(100*rem(tau,1))];

if ~exist([pwd '\' geomName '\' folderName], 'dir')
    mkdir([pwd '\' geomName '\' folderName]);
end

% cd([pwd '\' geomName '\' folderName]);
test.post_process([pwd '\' geomName '\' folderName]);

copyfile('woodward_test.m',[pwd '\' geomName '\' folderName '\case_setup.m']);

% cd('../../');