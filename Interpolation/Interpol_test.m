clear;
close all;

%% Setup

% test domain
dx = 0.5;
start_val = 0;
end_val = 2;

%% test inputs

coeffs = [0.1, -0.5, 6, -8, 0.7, 0.8];
res = 10;
f_in = @(x) polyval(coeffs, x);
fx_in = @(x) polyval(coeffs(1:end-1).*(length(coeffs(1:end-1)):-1:1), x);
fxx_in = @(x) polyval(coeffs(1:end-2).*(length(coeffs(1:end-1)):-1:2).*(length(coeffs(1:end-2)):-1:1), x);

xin = linspace(0, end_val, (end_val - start_val)/dx+1);
yin = f_in(xin);
k1 = fx_in(xin(1));%(cos(xin(1)))^2 - (sin(xin(1)))^2;
kn = fx_in(xin(end));%(cos(xin(end)))^2 - (sin(xin(end)))^2;
g1 = fxx_in(xin(1));%-4*sin(xin(1))*cos(xin(1));
gn = fxx_in(xin(end));%-4*sin(xin(end))*cos(xin(end));

SOR = 0.54;

%% Run quintic interpolation
[xout, yout, count, resid_k, resid_g, res1k, res1g] = quintic_interp(xin, yin, k1, kn, g1, gn, SOR, res);

figure();
scatter(xin, yin);
hold on;
plot(xout, yout, '*');
hold on;

y_test = f_in(xout);
res = max(abs(y_test - yout));
plot(xout, y_test);
legend('Control Points', 'Quintic Interpolations', 'Original Function');

title('Quintic Interpolation Verification');

saveas(gcf, 'quintic.png');

fprintf('Iteration Count: %i\nResidual (k): %0.5f\nResidual (g): %0.5f\n', count, max(resid_k), max(resid_g));
fprintf('First residual for k''s: %0.5f\n', max(res1k));
fprintf('First residual for g''s: %0.5f\n', max(res1g));
fprintf('Errors in Outputted Plots: %0.6f\n', res);


%% Bicubic Interpolation

% Define Grid
dx = 0.5;
dy = 0.25;
x_start = -1.0;
x_end = 2.0;
y_start = -1.0;
y_end = 2.0;

xvals = linspace(x_start, x_end, (x_end - x_start)/dx + 1);
yvals = linspace(y_start, y_end, (y_end - y_start)/dy + 1);

[XX, YY] = meshgrid(xvals, yvals);
[y_ct, x_ct] = size(XX);

res = 5; % resolution of the plots

% original function and analytical derivatives
fx = @(x, y) f_in(x) + f_in(y);
dfdx = @(x, y) fx_in(x);
dfdy = @(x,y) fx_in(y);
dfdxy = @(x,y) 0;

% Original solution and boundary conditions
ZZ = fx(XX, YY);%sin(XX).*cos(YY);
dfdx_1 = dfdx(XX(:,1),YY(:,1));
dfdx_n = dfdx(XX(:,end), YY(:,end));
dfdy_1 = dfdy(XX(1,:),YY(1,:));
dfdy_n = dfdy(XX(end,:), YY(end,:));
dfdxy_x1 = dfdxy(XX(:,1),YY(:,1));
dfdxy_xn = dfdxy(XX(:,end),YY(:,end));
dfdxy_y1 = dfdxy(XX(1,:),YY(1,:));
dfdxy_yn = dfdxy(XX(end,:),YY(end,:));

% get and plot results
[xout, yout, zout] = bicubic_interp(XX, YY, ZZ, dfdx_1, dfdx_n, dfdy_1, dfdy_n, dfdxy_x1, dfdxy_xn, dfdxy_y1, dfdxy_yn, res);

figure();

figure();
plot3(xout,yout,zout, '.');
title('Bicubic Spline Verification')
saveas(gcf, 'bicubic.png')
