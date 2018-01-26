clc;
clear;
close all;

%% Load Data

potential = load('potential_14.mat', 'OUT', 'GR');
euler = load('eulerIsen.mat','OUT', 'GR');

%% Plot Bow Shocks on X-axis

figure();
xMax = -3;
xRangePot = potential.GR.XX(:,end)>=xMax;
xRangeEul = euler.GR.XX(:,end)>=xMax;
plot(potential.GR.XX(xRangePot,end),potential.OUT.Uvals(xRangePot,end,1,end));
hold on;
plot(euler.GR.XX(xRangeEul,end), euler.OUT.Uvals(xRangeEul,end,1,end));
set(gca,'YDir','reverse');
title('Comparison of Bow Shocks - Density');
xlabel('X');
ylabel('\rho');
legend('Full Potential', 'Isentropic Euler', 'Location', 'Best');
set(gca, 'FontSize', 16);
saveas(gcf, 'density_shock', 'pdf');

%% Plot Contour Fields
figure();
contourf(potential.GR.XX, potential.GR.YY, potential.OUT.Uvals(:,:,1,end), 50);
axis equal;
title('Density Field - Potential Flow');
colorbar;
set(gca,'FontSize', 16);
saveas(gcf, 'density_field_pot', 'pdf');