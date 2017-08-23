clc;
close all;
clear;

case_name = 'test_object';

%% Define grid_lims

dx = 0.05;
dy = 0.08;

% Field Axis Values - body fitted grid
grid_lims = [   0,  15.0, dy; % y-range
                -5, 10,    dx];... 

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
CT.eps_s = 0.07525; % spatial diffusion term
CT.eps_t = 0.013;  % time diffusion term
CT.tol = 1e-6;
CT.dt = 0.1;
CT.iter_min = 300;
CT.CFL_on = 1;
CT.use_1visc = 1;
CT.is_polar = 0;
CT.case_name = 'cylinder_vec_1visc';

%% Fluid Params
FL.M0 = 1.01;
FL.gam = 1.4;

%% Test genBC

BC_setup.phys_types = {  {'inlet'},...
                {'far-field'},...
                {'outlet'},...
                {'wall'}};
BC_setup.vals = {{{1,1,0}},...
        {{1,1,0}},...
        {{0,0,0}},...
        {{0,0,dyBdx}}};
    
BC_setup.ranges = { {},...
                    {},...
                    {},...
                    {}};
        
test = eulerIsentropicField(CT, grid_lims(:,:,1), FL, BC_setup);
tic;
sim_time = 0;
while test.checkConvergence
    test = test.timeStep_tl;
    test = test.updateVals;
    
    sim_time(end+1) = test(1).CT.dt+sim_time(end);
    if abs(sim_time(end) - 2) < 1e-5
       pause('on'); 
    end
end

%% Post Process
folderName = ['M_' num2str(FL.M0)];
geomName = [case_name];
sizeName = [num2str(test(1).GR.d1) '_' num2str(test(1).GR.d2)];
dirName = [pwd '\' geomName '\' sizeName '\' folderName];


if ~exist(dirName, 'dir')
    mkdir(dirName);
end

% cd([pwd '\' geomName '\' folderName]);
test.post_process(dirName);

copyfile('woodward_test.m',[dirName '\case_setup.m']);