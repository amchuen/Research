close all;
clear;
clc;

%% Setup

% Inputs
nvals = [51, 101, 201, 401, 801];
u_0 = 1;
u_n = cosh(1);
dx = (1 - 0)./(nvals-1);

% Results storage
x_i = cell(1,length(nvals));
u_i = cell(1,length(nvals));

%% Solve and Plot 

figure();
error = zeros(size(nvals));
for i = 1:length(nvals)
    % Solve    
    x_i{i} = linspace(0,1, nvals(i));
    u_i{i} = ODEsol_Splines(x_i{i}, u_0, u_n);
    
    % Get Errors
    error(i) = trapz(x_i{i}, abs(u_i{i} - cosh(x_i{i})));
    
end

% Plot error relationship
plot(dx, error, '-+');
xlabel('Domain Spacing');
ylabel('Error');


