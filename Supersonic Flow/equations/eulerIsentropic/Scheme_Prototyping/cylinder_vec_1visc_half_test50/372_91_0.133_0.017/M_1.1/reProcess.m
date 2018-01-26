clc;
clear;
close all;

%% Load data files

load('results');

%% Round out numbesr

GR.r_cyl = 0.5;

%% Plot Results


figure();
contourf(GR.XX, GR.YY, round(EE.fv(:,:,1,end),3), 5);axis equal;
title('Density, \rho');
xlabel('X');
ylabel('Y');
colorbar;
saveas(gcf, [ 'density'],'pdf');
saveas(gcf, [ 'density']);

cpVals = 1-q2_ij;
figure();contourf(GR.XX, GR.YY, round(cpVals, 3), 5);
axis equal;
title(['Coefficient of Pressure, M0 =' num2str(M0)]);
colorbar;
saveas(gcf, [ 'cp_contour'],'pdf');
saveas(gcf, [ 'cp_contour']);
