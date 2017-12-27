clc;
clear;
close all;

%% Initialize

dt = 0.01;
time = 0:dt:100;
iterateFlag = 1;

%% Generate Equations
%     qFunc = @(vMat, t)  [   -funcX(vMat(2,1),vMat(2,2), t), -funcY(vMat(2,1),vMat(2,2), t)];

delta = 0.22;
omega = 1;
beta = 1;
alpha = 1;
fVal = 0.3;

qFunc = @(vMat, t)  alpha.*vMat.^2 - beta;
pFunc = @(vMat, t)  delta;
rFunc = @(vMat, t)  fVal.*cos(omega .*t);

%% Solve using cubic2

figure(1);
figure(2);
figure(3);

v0 = [1;0]; % [xdot, ydot; x, y]
vCubic = cubicTimeStep2(time, pFunc, qFunc, rFunc, v0);

%% Solve using Explicit RK4

v_duff = @(vVec, t) fVal.*cos(omega.*t) - delta*vVec(1) + beta.*vVec(2) - alpha.*vVec(2)^3;
vFunc_duff = @(vVec, t) [   v_duff(vVec, t);...
                    vVec(1)];
       
vRK4_duff = rk4(time, vFunc_duff, v0);

%% Solve using cubic1

vCubic_duff = cubicTimeStep(time, v_duff, v0, 1);

%% Compare Results

figure(1);
hold on;
plot(time, reshape(vCubic(2,1,:), 1, length(vCubic)));
plot(time, vRK4_duff(2,:), '--');
plot(time, vCubic_duff(2,:), '-.');