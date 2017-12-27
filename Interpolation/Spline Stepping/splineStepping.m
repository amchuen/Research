clc;
clear;
close all;

%% Initialize

dt = 0.01;
time = 0:dt:200;
iterateFlag = 1;

% van der Pol equation coefficients
epsil = 10;
vQuad_duff = [1; 0]; %[xDot, x]
vQuad_vdp = vQuad_duff;
v_vdp = @(vVec, t) -epsil.*(vVec(2,:)^2 - 1).*vVec(1,:) - vVec(2,:);

% Duffing's Equation Coefficients
delta = 0.22;
omega = 1;
alpha = 1;
beta = 1;
ff = 0.3;

v_duff = @(vVec, t) ff.*cos(omega.*t) - delta*vVec(1) + beta.*vVec(2) - alpha.*vVec(2)^3;

%% Time-Step Using Quadratic Spline

for i = 1:length(time)-1
    % Calculate next step, generate quadratic spline at each time-step
    vQuad_duff(:,end+1) = [vQuad_duff(1,end)+v_duff(vQuad_duff(:,end), time(i)).*dt;...
                    vQuad_duff(2,end) + vQuad_duff(1,end).*dt + 0.5.*v_duff(vQuad_duff(:,end), time(i)).*dt^2];
    
    vQuad_vdp(:,end+1) = [vQuad_vdp(1,end)+v_vdp(vQuad_vdp(:,end), time(i)).*dt;...
                    vQuad_vdp(2,end) + vQuad_vdp(1,end).*dt + 0.5.*v_vdp(vQuad_vdp(:,end), time(i)).*dt^2];
end

%% TIme-Stepping Using Cubic Spline

vCubic_duff = cubicTimeStep(time, v_duff, vQuad_duff(:,1), 1);
vCubic_vdp = cubicTimeStep(time, v_vdp, vQuad_vdp(:,1), 1);

%% Time Stepping Using Quartic Spline

vQuart_duff = quarticTimeStep_iterQuart(time, v_duff, vQuad_duff(:,1), 100);
vQuart_vdp = quarticTimeStep_iterQuart(time, v_vdp, vQuad_vdp(:,1), 100);

%% Time-Step with RK4

vFunc_duff = @(vVec, t) [   v_duff(vVec, t);...
                    vVec(1)];
vFunc_vdp =  @(vVec, t) [   v_vdp(vVec, t);...
                    vVec(1)];
        
vRK4_duff = rk4(time, vFunc_duff, vQuad_duff(:,1));
vRK4_vdp = rk4(time, vFunc_vdp, vQuad_duff(:,1));

%% Compare Results

figure('pos', [10 10 900 600]);
subplot(3,1,1);
plot(time, vQuad_duff(2,:));
hold on;
plot(time, vRK4_duff(2,:), '-.');
title('Quadratic Spline Time-stepping');
legend('Quad', 'RK4');
% xlabel('t');
ylabel('X');
subplot(3,1,2);
plot(time, vCubic_duff(2,:)); hold on;plot(time, vRK4_duff(2,:), '-.');
title('Cubic Spline Time-stepping');
legend('Quartic', 'RK4');
% xlabel('t');
ylabel('X');
subplot(3,1,3);
plot(time, vQuart_duff(2,:)); hold on;plot(time, vRK4_duff(2,:), '-.');
title('Quartic Spline Time-stepping');
legend('Quartic', 'RK4');
xlabel('t');
ylabel('X');
suptitle('Duffing Equation');
% saveas(gcf, 'XvsT', 'pdf');

figure('pos', [10 10 900 600]);
subplot(3,1,1);
plot(time, vQuad_vdp(2,:));
hold on;
plot(time, vRK4_vdp(2,:), '-.');
title('Quadratic Spline Time-stepping');
legend('Quad', 'RK4');
% xlabel('t');
ylabel('X');
subplot(3,1,2);
plot(time, vCubic_vdp(2,:)); hold on;plot(time, vRK4_vdp(2,:), '-.');
title('Cubic Spline Time-stepping');
legend('Quartic', 'RK4');
% xlabel('t');
ylabel('X');
subplot(3,1,3);
plot(time, vQuart_vdp(2,:)); hold on;plot(time, vRK4_vdp(2,:), '-.');
title('Quartic Spline Time-stepping');
legend('Quartic', 'RK4');
xlabel('t');
ylabel('X');
suptitle('van der Pol Equation');
% saveas(gcf, 'XvsT', 'pdf');

figure();plot(vQuart_vdp(2,:), vQuart_vdp(1,:));
hold on;
plot(vRK4_vdp(2,:), vRK4_vdp(1,:), '--');
title('van der Pol (Quartic)');
legend('Quad', 'RK4');
xlabel('X');
ylabel('$\dot{X}$', 'Interpreter', 'latex');
% saveas(gcf, 'XDotvsX_quad', 'pdf');

figure();plot(vQuart_duff(2,:), vQuart_duff(1,:));
hold on;
plot(vRK4_duff(2,:), vRK4_duff(1,:), '--');
title('Duffing (Quartic)');
legend('Quartic', 'RK4');
xlabel('X');
ylabel('$\dot{X}$', 'Interpreter', 'latex');
% saveas(gcf, 'XDotvsX_quartic', 'pdf');

%% Error Checking

fprintf('Duffing Equation:\n');
fprintf('Error of Quad Spline: %0.5e\n', trapz(time, abs(vQuad_duff(2,:)-vRK4_duff(2,:))));
fprintf('Error of Cubic Spline: %0.5e\n', trapz(time, abs(vCubic_duff(2,:)-vRK4_duff(2,:))));
fprintf('Error of Quartic Spline: %0.5e\n', trapz(time, abs(vQuart_duff(2,:)-vRK4_duff(2,:))));

fprintf('\nvan der Pol Equation:\n');
fprintf('Error of Quad Spline: %0.5e\n', trapz(time, abs(vQuad_vdp(2,:)-vRK4_vdp(2,:))));
fprintf('Error of Cubic Spline: %0.5e\n', trapz(time, abs(vCubic_vdp(2,:)-vRK4_vdp(2,:))));
fprintf('Error of Quartic Spline: %0.5e\n', trapz(time, abs(vQuart_vdp(2,:)-vRK4_vdp(2,:))));