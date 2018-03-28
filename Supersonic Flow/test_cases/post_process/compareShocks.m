clc;
clear;
close all;

%% Load Data

potential14 = load('potential_14.mat', 'OUT', 'GR', 'FL');
potential11 = load('potential_11.mat', 'OUT', 'GR', 'FL');
euler14 = load('eulerIsen_14.mat','OUT', 'GR','FL');
euler11 = load('eulerIsen_11.mat','EE', 'GR', 'M0', 'gam');

%% Plot Bow Shocks on X-axis

figure();
xMax = -13;
xRangePot = potential14.GR.XX(:,end)>=xMax;
xRangePot11 = potential11.GR.XX(:,end)>=xMax;
xRangeEul = euler14.GR.XX(:,end)>=xMax;
xRangeEul11 = euler11.GR.XX(:,end)>=xMax;
plot(potential14.GR.XX(xRangePot,end),potential14.OUT.Uvals(xRangePot,end,1,end));
hold on;
plot(euler14.GR.XX(xRangeEul,end), euler14.OUT.Uvals(xRangeEul,end,1,end));
hold on;
plot(euler11.GR.XX(xRangeEul11,end), euler11.EE.fv(xRangeEul11,end,1,end));
hold on;
plot(potential11.GR.XX(xRangePot11,end),potential11.OUT.Uvals(xRangePot11,end,1,end));
set(gca,'YDir','reverse');
title(['Comparison of Bow Shocks - Density ' num2str(potential14.FL.M0)]);
xlabel('X');
ylabel('\rho');
legend('Full Potential (M0 = 1.4)', 'Isentropic Euler (M0 = 1.4)', 'Isentropic Euler (M0=1.1)', 'Location', 'Best');
set(gca, 'FontSize', 16);
saveas(gcf, 'density_shock', 'pdf');

%% Plot Contour Fields
figure();
contourf(potential14.GR.XX, potential14.GR.YY, potential14.OUT.Uvals(:,:,1,end), 10);
axis equal;
title('Density Field - Potential Flow');
colorbar;
set(gca,'FontSize', 16);
saveas(gcf, 'density_field_pot', 'pdf');