clc;
clear;
close all;

%% Initialize Grid and Field Values

xx = linspace(-1,1, 201);
dx = xx(2) - xx(1);
dt = 0.4*dx;
dtLast = dt;
cfl = 1;

UU = zeros(size(xx));
UU(xx <= 0) = 1;
UU(xx > 0) = -1;
UU = repmat(UU,3,1);

%% Important Functions

% visc_x = 0.2;
% d2U = @(U) (U(3:end) - 2.*U(2:end-1) + U(1:end-2));%./(dx^2);
% flux = @(U) U(2:end-1) .* (U(3:end) - U(1:end-2))./(2*dx)+...
%             (visc_x.*([d2U(U(2:end)),0] - d2U(U)) - visc_x.*(d2U(U) - [0,d2U(U(1:end-1))]))./dx - ...
%             (VRvisc(U(2:end),dx).*(U(3:end) - U(2:end-1)) - VRvisc(U(1:end-1),dx).*(U(2:end-1)-U(1:end-2)))./(dx^2);
        
%% Time Step

tol = 1e-6;
flux = fluxFunc(UU(end,:),dx);
res = norm(max(abs(flux)));

figure(1);
resPlot = semilogy(res);
title('Residual');
movegui(gcf, 'west');

figure(2);
Ux = plot(xx, UU(3,:));
movegui(gcf, 'east');
ratio = 8;

while res(end) > tol
   
    % Update Values
    UU(1:2,:) = UU(2:3,:);
    
    % Correct Time Steps
    if (max(abs(UU(2,:)))*dt/dx ~= cfl) %|| (dt > (dx)/(8*visc_x))
%         dt = min(cfl.*dx/max(abs(UU(2,:))), dx/(8*visc_x));
        dt = (cfl.*dx/max(abs(UU(2,:))));
        if abs(log10(dt/dtLast)) > 0.301 %1e-4
            dtLast = dt;
            fprintf('Time-step changing!\nNew time step: %0.5e\n', dt);
        end
    end
    
    % Calculate next time step
%     UU(3,2:end-1) = UU(2,2:end-1) - dt.*fluxFunc(UU(2,:),dx);
%     beta = 4.*visc_x*(dt^2)/(ratio*dx^2);
    beta = 0.5.*(VRvisc(UU(3,2:end),dx)+VRvisc(UU(3,1:end-1),dx)).*(dt/dx)^2;
    UU(3,2:end-1) = ((2.*beta./(dt^2)).*UU(2,2:end-1)-(beta./(dt^2)-0.5./dt).*UU(1,2:end-1)-flux)./(beta./(dt^2)+0.5/dt);
    
    % Update Residuals
    flux = fluxFunc(UU(3,:),dx);
    res(end+1) = norm(max(abs(flux)));
    resPlot.YData(end+1) = res(end);
    Ux.YData = UU(3,:);
    drawnow;
    
end

%% Post Process

figure(1);
saveas(gcf, 'residual.pdf');

figure();plot(xx, UU(3,:), 'o-');
saveas(gcf, '1D_1Eqn_results', 'pdf');