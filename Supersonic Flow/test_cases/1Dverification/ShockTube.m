clc;
clear;
close all;

%% Generate Grid

% dx = 0.0025;
% xx = (0:dx:1)';
xx = linspace(0,1,1001);
dx = xx(2) - xx(1);
dt = 0.1;

%% Fluid Properties
gam = 1.4;
cfl = 0.95;

flux =@(U) [ U(2,:,2);...
            0.5.*(3-gam).*(U(2,:,2).^2) ./ U(1,:,2) + (gam-1).*U(3,:,2);...
            -0.5.*(gam-1).*(U(2,:,2).^3)./(U(1,:,2).^2) + gam.*U(2,:,2).*U(3,:,2)./U(1,:,2)];

%% Initialize Field Values

UU = ones([3, length(xx)]);

% density
UU(1,:) = (xx<=0.5)+0.125.*(xx>0.5); 

% velocity
UU(2,:) = 0; 

%Energy
UU(3,:) = 2.5.*(xx<=0.5)+0.25.*(xx>0.5); 

UU = repmat(UU,1,1,3);

%% Setup Time Simulation
time = 0;
tEnd = 10*5;

%% Run Simulation

while length(time) < 500
    
    % Update Field Values
    UU(:,:,1:2) = UU(:,:,2:3);
    
    % Update Flux
%     PP = (gam - 1).*(UU(:,3,2) - 0.5.*UU(:,2,2).^2./UU(:,1,2));
    FF =  flux(UU);
    
    CFL1 = dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    CFL1 = max(abs(CFL1(:)));
    CFL2 = 2.*(UU(2,2:end-1,end)./UU(1,2:end-1,end)).*dt./(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    CFL2 = max(abs(CFL2(:)));
    cflFactor = 1;
    if (~isempty(CFL1) && ((CFL1 > cfl)||(CFL2 ~=cfl))) && any(UU(2,:,2) > 0)
        cflFactor = min(cfl/CFL1,cfl/CFL2);
        dt = dt.*cflFactor;
        
        if abs(log10(cflFactor)) > 0.1
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps\n');
            fprintf('New time step:%0.5e\n', dt);
        end
        
    end    
    
    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (UU(:,2:end-1,1).*(1 - dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx))) - dt.*(FF(:,3:end) - FF(:,1:end-2))./dx + 2.*dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx).*UU(:,3:end,2) + epsFunc(UU(:,1:end-1,2), dx).*UU(:,1:end-2,2)))./(1 + dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx)));
    time(end+1) = time(end)+dt;
    
end

%% Post Process

figure(1);
for i = 1:size(UU,1)
    if i == 1
        plot(xx, UU(i,:,end));
    else
        plot(xx, UU(i,:,end)./UU(1,:,end));
    end
    hold on;
end

PP = (gam - 1).*(UU(3,:,3) - 0.5.*UU(2,:,3).^2./UU(1,:,3));
plot(xx,PP);

hold off;
title('Sod Shock Tube');
legend('\rho', 'u', 'e', 'p', 'Location', 'BestOutside');
set(gca,'FontSize', 16);
saveas(gcf, 'sodShockTube.pdf');