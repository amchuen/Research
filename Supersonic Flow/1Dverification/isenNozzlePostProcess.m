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
plot(gam_list, xShock, 'o');

pCoeffs = polyfit(gam_list, xShock, 3);
hold on;
x1 = linspace(gam_list(1), gam_list(end));
y1 = polyval(pCoeffs, x1);
plot(x1, y1);
title('\gamma vs Shock Location');
xlabel('\gamma');
ylabel('X');

txt1 = ['y =' num2str(pCoeffs(end)) 'x^3' num2str(pCoeffs(end-1)) 'x^2 +' num2str(pCoeffs(end-2)) 'x' num2str(pCoeffs(1))];
text(1.6, 0.8, txt1, 'FontSize', 13);
set(gca, 'FontSize', 14);
saveas(gcf, [dirName '\gamma_curveFit'], 'pdf');
saveas(gcf, [dirName '\gamma_curveFit'], 'fig');

%% Plot Shocks

figure();
for i = 1:7
    plot(xx, results(i).UU(2,:,3)./results(i).UU(1,:,3), '*-');
    hold on;
end
plot(xShock, uShock, 'o');

save([dirName '\' 'postProcess'], 'xShock', 'uShock', 'gam_list');