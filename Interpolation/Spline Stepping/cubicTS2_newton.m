clc;
clear;
close all;

%% Initialize

dt = 0.01;
time = 0:dt:50;
iterateFlag = 1;

mu = [0.4, 0.5, 1, 0.6];

figure(1);
figure(2);
figure(3);

for i = 1:length(mu)
    funcX = @(x, y, t) mu(i)./(x^2 + y^2).^(1.5); % <- input w^2, not q(x)*y!
    funcY = @(x, y, t) mu(i)./(x^2 + y^2).^(1.5);

%     qFunc = @(vMat, t)  [   -funcX(vMat(2,1),vMat(2,2), t), -funcY(vMat(2,1),vMat(2,2), t)];
    qFunc = @(vMat, t)  [   funcX(vMat(1),vMat(2), t), funcY(vMat(1),vMat(2), t)];
    pFunc = @(vMat, t)  [   0, 0];
    rFunc = pFunc;
    %                         vVec(1,1),  vVec(2,1)];
    v0 = [0, 1; 1, 0]; % [xdot, ydot; x, y]
    vCubic = cubicTimeStep2(time, pFunc, qFunc, rFunc, v0);
    figure(1);
    subplot(length(mu),1,i);
    hold on;
    plot(reshape(vCubic(2,1,:), 1, length(vCubic)), reshape(vCubic(2,2,:), 1, length(vCubic)));
    
    rk4Func = @(vVec, t) [  -funcX(vVec(2),vVec(4), t).*vVec(2);...
                            vVec(1);
                            -funcY(vVec(2),vVec(4), t).*vVec(4);...
                            vVec(3)];
    vRK4 = rk4(time, rk4Func, [0;1;1;0]);

    hold on;
    plot(vRK4(2,:), vRK4(4,:), '--');
    title(['Position for mu = ', num2str(mu(i))]);
    
    Hcubic = reshape(0.5.*sum(vCubic(1,:,:).^2,2) - 1./sqrt(sum(vCubic(2,:,:).^2, 2)), 1, length(time));
    Mcubic = reshape(vCubic(2,1,:).*vCubic(1,2,:) - vCubic(2,2,:).*vCubic(1,1,:),1, length(time));
    Hrk4 = 0.5.*sum(vRK4(1:2:3,:).^2,1) - 1./sqrt(sum(vRK4(2:2:4,:).^2, 1));
    MRK4 = prod(vRK4(2:3,:),1) - prod(vRK4(1:3:4,:),1);
    figure(2);
    subplot(length(mu),1,i);
    hold on;
    plot(time, Hcubic); hold on; plot(time, Hrk4, '--');
    title(['mu = ', num2str(mu(i))]);
%     axis equal;
    suptitle('Energy vs Time');
    
    figure(3);
    subplot(length(mu),1,i);
    hold on;
    plot(time, Mcubic); hold on; plot(time, MRK4, '--');
    title(['mu = ', num2str(mu(i))]);
%     axis equal;
    suptitle('Moment vs Time');
    
    fprintf('Digits of Accuracy Lost: %0.2f\n', abs(floor(log10(eps))) - abs(floor(log10(max(abs(vRK4(2,:) - reshape(vCubic(2,1,:), 1, length(vCubic))))))));
end