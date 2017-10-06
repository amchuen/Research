clc;
close all;
clear;

case_name = 'test_object';

%% Define grid_lims

dx = 0.01;
dy = 0.01;

% Field Axis Values - body fitted grid
grid_lims = [   0+(0.5*dy),  1.0 - (0.5*dy), dy; % y-range
                0+(0.5*dx), 0.6 - (0.5*dx),    dx];
                
grid_lims(:,:,2) = [   0.2+(0.5*dy),  1.0-(0.5*dy), dy; % y-range
                        0.6-(0.5*dx),  3.0-(0.5*dx),  dx];            

%% CTRL
CT.eps_s = 0.002525; % spatial diffusion term
CT.eps_t = 0.013;  % time diffusion term
CT.tol = 1e-6;
CT.dt = 0.006;
CT.iter_min = 50000;
CT.CFL_on = 0;
CT.maxCFL = 1.2;
CT.use_1visc = 1;
CT.is_polar = 0;
CT.case_name = 'cylinder_vec_1visc';

%% Fluid Params
FL.M0 = 3.0;
FL.gam = 1.4;

%% Test genBC

BC_setup.phys_types = {  {'inlet'},...
                {'wall'},...
                {'wall', 'patch'},...
                {'wall'}};
BC_setup.vals = {{{1,1,0}},...
        {{0,0,0}},...
        {{0,0,0},{1,1,0}},...
        {{0,0,0}}};
    
BC_setup.ranges = { {},...
                    {},...
                    {[0, 0.2], [0.2, 1]},...
                    {}};
BC_setup(2).phys_types = {  {'patch'},...
                {'wall'},...
                {'outlet'},...
                {'wall'}};
BC_setup(2).vals = {{{1,1,0}},...
        {{0,0,0}},...
        {{0,0,0}},...
        {{0,0,0}}};
    
BC_setup(2).ranges = { {[0.2, 1]},...
                    {},...
                    {},...
                    {}};                
        
test = eulerIsentropicField(CT, grid_lims(:,:,1), FL, BC_setup(1));
test(2) = eulerIsentropicField(CT, grid_lims(:,:,2), FL, BC_setup(2));
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