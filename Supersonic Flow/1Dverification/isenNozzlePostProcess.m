clc;
clear;
close all;

dirName = 'isentropicNozzle';

%% Load Data

xx = linspace(0.5,1,151);
dx = xx(2) - xx(1);

for i = 1:7
    results(i) = load([dirName '\' 'gam_' num2str(i)]);
end

%% Calculate Shock Location

gam_list = [];
xShock = [];

for i = 1:7
   results(i).dudx = abs(diff(results(i).UU(2,:,3)./results(i).UU(1,:,3),1,2));
   results(i).indShock = find(results(i).dudx == max(results(i).dudx));
   xShock(i) = 0.5*(xx(results(i).indShock) + xx(results(i).indShock+1));
   uShock(i) = 0.5*(results(i).UU(2,results(i).indShock,3)./results(i).UU(1,results(i).indShock,3) + results(i).UU(2,results(i).indShock+1,3)./results(i).UU(1,results(i).indShock+1,3));
   gam_list(i) = results(i).gam;
end

%% Plot Data and Curve-Fit!
figure();
plot(gam_list, xShock, 'o-');

%% Plot Shocks

figure();
for i = 1:7
    plot(xx, results(i).UU(2,:,3)./results(i).UU(1,:,3));
    hold on;
end
plot(xShock, uShock, 'o');