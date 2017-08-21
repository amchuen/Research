clc;
close all;
clear;

%% Define grid_lims

dx = 0.01;
dy = 0.01;

% Field Axis Values - body fitted grid
grid_lims = [   0+0.5*dy,  0.2-0.5*dy, dy; % y-range
                0.0+(0.5*dx), 0.6-(0.5.*dx),    dx];... 

grid_lims(:,:,2) = [    0.2+0.5.*dy-mod(0.2,dy),  1.0, dy; % y-range
                        0.0+(0.5*dx), 0.6-(0.5.*dx),    dx];... 
                        
grid_lims(:,:,3) = [    0.2+0.5.*dy-mod(0.2,dy),  1.0, dy;
                        0.6+(0.5.*dx)-mod(0.6,dx), 3.0, dx];

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

%% CTRL
CT.eps_s = 0.001;% 0.07525; % spatial diffusion term
CT.eps_t = 0.013;  % time diffusion term
CT.tol = 1e-7;
CT.dt = 0.1;
CT.iter_min = 300;
CT.CFL_on = 1;
CT.use_1visc = 1;
CT.is_polar = 0;
CT.case_name = 'cylinder_vec_1visc';

%% Fluid Params
FL.M0 = 3.0;
FL.gam = 1.4;

%% Test genBC

BC_setup.phys_types = {  {'inlet'},...
                {'patch'},...
                {'wall', 'outlet'},...
                {'wall'}};
BC_setup.vals = {{{1,1,0}},...
        {{1,1,0}},...
        {{0,0,0}},...
        {{0,0,dyBdx}}};
    
BC_setup.ranges = {   {},...
            {},...
            {[0, 0.1], [0.1, 0.2]},...
            {}};
        
test = eulerIsentropicField(CT, grid_lims(:,:,1), FL, BC_setup);
% BC_setup = {phys_types, vals, range};

% BC = genBC(phys_types, vals, range);