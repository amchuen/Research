clc;
clear;
close all;

%% Euler Equations of Rigid Body Rotation (with Free Torque and Beta = 0)
dt = 0.01;
time = linspace(0, 80, 5001);
% time = 0:dt:80;
iterateFlag = 1;

I1 = 3; I2 = 2; I3 = 1;
beta = 0.1;
funcW1 = @(wVec, w_t, t)    ((wVec(3)^2)*(I2-I3).*(I3-I1)./(I1*I2) +...
                            (wVec(2)^2)*(I2-I3).*(I1-I2)./(I1*I3)).*wVec(1) - ...
                            beta .* w_t(1);
                        
funcW2 = @(wVec, w_t, t)    ((wVec(1)^2)*(I3-I1).*(I1-I2)./(I2*I3) +...
                            (wVec(3)^2)*(I2-I3).*(I3-I1)./(I1*I2)).*wVec(2) - ...
                            beta .* w_t(2);
                        
funcW3 = @(wVec, w_t, t)    ((wVec(2)^2)*(I1-I2).*(I2-I3)./(I1*I3) + ...
                            (wVec(1)^2)*(I1-I2).*(I3-I1)./(I2*I3)).*wVec(3) - ...
                            beta .* w_t(3);

% Cubic Spline Input
dFunc_rot = @(wMat, t) [funcW1(wMat(2,:),wMat(1,:),t), funcW2(wMat(2,:),wMat(1,:),t), funcW3(wMat(2,:),wMat(1,:),t)];

% RK4 Input
rotRK4 = @(wVec, t) [   funcW1(wVec(4:6), wVec(1:3), t);...
                        funcW2(wVec(4:6), wVec(1:3), t);
                        funcW3(wVec(4:6), wVec(1:3), t);
                        wVec(1:3)];
                    
w0 = [  0, 0, 0;...
        1, 0.01, 0.01];

% international cfd conference

wCubic = cubicTimeStep(time, dFunc_rot, w0, 1);
wRK4 = rk4(time, rotRK4, [0;0;0;1;0.01;0.01]);

%% Plot Angular Motion
motionCubic = reshape(wCubic(2,:,:), size(wCubic,2), size(wCubic,3));
dmotCubic = reshape(wCubic(1,:,:), size(wCubic,2), size(wCubic,3));
ddmotCubic = zeros(size(dmotCubic));
for i = 1:length(time)
   ddmotCubic(:,i) = dFunc_rot(wCubic(:,:,i),time(i))';    
end
marktype = {'-', '--', 'k-.', '-.'};
figure(1);
suptitle('Time History of Angular Velocity, \Deltat = 0.016');
figure(2);
figure(3);
for i = 1:size(motionCubic,1)
    figure(1);
    ax(i) = subplot(size(motionCubic, 1), 1, i);
    plot(time, motionCubic(i,:), marktype{i}); hold on;
    ylabel(['\omega_' num2str(i)], 'FontSize', 14);
%     if i < size(motionCubic,1)
%         set(gca,'XTick',[]);
%     end
    figure(2);
    plot(time, dmotCubic(i,:), marktype{i}); hold on;
    figure(3);
    plot(time, ddmotCubic(i,:), marktype{i}); hold on;
end
figure(1);
xlabel('Time','FontSize', 14);
linkaxes(ax,'x');
saveas(gcf,'eulerSoln', 'pdf');

%% Plot Kinetic Energy
energyCubic = (0.5.*reshape(wCubic(2,:,:).^2, size(wCubic,2), size(wCubic,3))'*[I1;I2;I3])';
energyRK4 = 0.5.*((wRK4(4:end,:).^2)'*[I1;I2;I3])';

figure();plot(time, energyCubic, 'LineWidth', 0.75);
hold on;
plot(time, energyRK4, '*-.', 'MarkerIndices', 1:50:length(time));
title('Kinetic Energy');
xlabel('time');
ylabel('Energy');
legend('Cubic Interpolation', 'RK4');
saveas(gcf,'eulerEnergy', 'pdf');

errEnergy = trapz(time, abs(energyCubic - energyRK4));
test = simpInteg(time, abs(energyCubic - energyRK4));

figure();semilogy(time, abs(energyCubic - energyRK4));
title({'Time History of Energy Error','(Rigid Body Rotation, $\Delta t = 0.016$)'}, 'FontSize', 18, 'Interpreter', 'Latex');
xlabel('Time', 'FontSize', 14);
ylabel('Error (Energy)', 'FontSize', 14);
saveas(gcf, 'eulerEnergyError','pdf');

%% Plot Angular Momentum

angCubic = (reshape(wCubic(2,:,:).^2, size(wCubic,2), size(wCubic,3))'*([I1;I2;I3].^2))';
angRK4 = ((wRK4(4:end,:).^2)'*([I1;I2;I3].^2))';

figure();plot(time, angCubic);
hold on;
plot(time, angRK4, '--');
title('Angular Momentum');