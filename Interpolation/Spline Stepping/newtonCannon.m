clc;
clear;
close all;

%% Initialize

kdt = [2, 1, 0.5];
dt = 0.005;

iterateFlag = 1;

mu = [0.4, 0.5, 1, 0.6];

figure(1);
figure(2);
figure(3);
figure(4);

epsilon = [mu(1), 0, -mu(3)./(2*1), -mu(end)./(2.*3)];

Hconvergence = zeros(length(mu), length(kdt));

for ii = 1:length(kdt)
    
    time = 0:kdt(ii)*dt:50;
    timeRK4 = 0:kdt(ii)*dt:50;
    fprintf('\nStep Size:%0.4f\n', kdt(ii)*dt);
    
    for i = 1:length(mu)
        funcX = @(x, y, t) -mu(i)*x./(x^2 + y^2).^(1.5);
        funcY = @(x, y, t) -mu(i)*y./(x^2 + y^2).^(1.5);

        dFunc = @(vVec, t)  [   funcX(vVec(2,1),vVec(2,2), t), funcY(vVec(2,1),vVec(2,2), t)];...
        %                         vVec(1,1),  vVec(2,1)];
        v0 = [0, 1; 1, 0]; % [xdot, ydot; x, y]
        vCubic = cubicTimeStep(time, dFunc, v0,1);

        rk4Func = @(vVec, t) [  funcX(vVec(2),vVec(4), t);...
                                vVec(1);
                                funcY(vVec(2),vVec(4), t);...
                                vVec(3)];
        vRK4 = rk4(timeRK4, rk4Func, [0;1;1;0]);
        
        % Calculate energy and momentum
        Hcubic = reshape(0.5.*sum(vCubic(1,:,:).^2,2) - 1./sqrt(sum(vCubic(2,:,:).^2, 2)), 1, length(time));
        Mcubic = reshape(vCubic(2,1,:).*vCubic(1,2,:) - vCubic(2,2,:).*vCubic(1,1,:),1, length(time));
        Hrk4 = 0.5.*sum(vRK4(1:2:3,:).^2,1) - 1./sqrt(sum(vRK4(2:2:4,:).^2, 1));
        MRK4 = prod(vRK4(2:3,:),1) - prod(vRK4(1:3:4,:),1);

        if kdt(ii) == 2
            figure(1);
%             subplot(length(mu),1,i);
            hold on;
            plot(reshape(vCubic(2,1,:), 1, length(vCubic)), reshape(vCubic(2,2,:), 1, length(vCubic)), 'LineWidth', 1.75);
%             hold on;
%             plot(vRK4(2,:), vRK4(4,:), '--');
%             title(['Position for mu = ', num2str(mu(i))], 'FontSize', 14);
            
            figure(2);
            subplot(length(mu),1,i);
            hold on;
            plot(time, Hcubic); hold on; plot(timeRK4, Hrk4, '--');
            title(['mu = ', num2str(mu(i))]);
            suptitle('Energy vs Time');
            
            figure(3);
            hold on;
            plot(time, abs(Hcubic - Hrk4));
            title('Energy Difference vs. Time -- Newton''s Cannonball');

            figure(4);
            subplot(length(mu),1,i);
            hold on;
            plot(time, Mcubic); hold on; plot(timeRK4, MRK4, '--');
            title(['mu = ', num2str(mu(i))]);
            suptitle('Moment vs Time');
        end
        
        Hconvergence(i,ii) = norm(Hcubic - Hrk4, inf);


        fprintf('Digits of Accuracy Lost: %0.2f\n', abs(floor(log10(eps)) - floor(log10(max(abs(Mcubic - 1))))));
    end
    
end

%% Post Process?

figure(1);
xlabel('X', 'FontSize', 14);
ylabel('Y', 'FontSize', 14);
legend({'$\mu = 0.4$', '$\mu = 0.5$', '$\mu = 1$', '$\mu = 0.6$'}, 'Interpreter', 'latex', 'Location', 'Best', 'FontSize', 16);
title({'Orbital Trajectories For Various', 'Gravity Parameters $\mu$, $\Delta t$ = 0.01'}, 'Interpreter', 'Latex', 'FontSize', 18);
saveas(gcf, 'newtonOrbit', 'pdf');


figure(2);
saveas(gcf, 'newtonEnergy','pdf');

figure(3);
xlabel('Time');
ylabel('Error (Energy)', 'FontSize', 14);
title({'Time History of Energy Error','(Rigid Body Rotation, $\Delta t = 0.01$)'}, 'FontSize', 18, 'Interpreter', 'Latex');
legend({'$\mu = 0.4$', '$\mu = 0.5$', '$\mu = 1$', '$\mu = 0.6$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'Best');
set(gca,'yscale', 'log');
saveas(gcf, 'newtonEnergyError', 'pdf');

figure();
for i = 1:length(mu)
   loglog(kdt.*dt, Hconvergence(i,:), 'LineWidth', 1.5);
   hold on;
end
set(gca,'XDir', 'reverse');
xlabel('Time Step-Size', 'FontSize', 14);
ylabel('Energy Error', 'FontSize', 14);
legend({'$\mu = 0.4$', '$\mu = 0.5$', '$\mu = 1$', '$\mu = 0.6$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'Best');
title('Energy Convergence (Newton''s Cannonball)', 'FontSize', 18);
saveas(gcf, 'newtonConvergence', 'pdf');