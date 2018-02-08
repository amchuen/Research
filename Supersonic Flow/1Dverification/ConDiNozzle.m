clc;
clear;
close all;

%% Generate Grid

% dx = 0.0025;
% xx = (0:dx:1)';
xx = linspace(0,1,101);
dx = xx(2) - xx(1);
dt = 0.1;

%% Nozzle Geometry
ind1 = xx>=0 & xx <=0.5;
ind2 = xx>0.5;
alpha = 0.4;
g_x = [9/4.*(2.*xx(ind1) - 1).^4 - 1.5.*(2.*xx(ind1)-1).^2 + 5/4, (1-alpha).*(9/4.*(2.*xx(ind2)-1).^4 - 1.5.*(2.*xx(ind2)-1).^2) + 5/4];
dgdx = [0, (g_x(3:end) - g_x(1:end-2))./(2.*dx), 0];

%% Fluid Properties
gam = 1.4;
cfl = 1;

flux =@(U) [ U(2,:,2);...
            0.5.*(3-gam).*(U(2,:,2).^2)./U(1,:,2) + (gam-1).*U(3,:,2);...
            -0.5.*(gam-1).*(U(2,:,2).^3)./(U(1,:,2).^2) + gam.*U(2,:,2).*U(3,:,2)./U(1,:,2)] .* g_x;

%% Initialize Field Values

% Inlet Conditions -> s0 = 0 (isentropic relations)
rho0 = (0.5.*(gam+1))^(1/(gam-1));
H0 = 3;
% u0 = sqrt(2*(H0 - (rho0^(gam-1))/(gam-1)));
u0 = 0;
% rho0 = ((gam-1).*(H0 - 0.5*u0^2)).^(1/(gam-1)); 
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;

% Exit Conditions
p_i = 0.9;
rho_i = 2*gam/(gam-1)*p_i;
E_i = p_i/((gam-1)*rho_i) + 0.5.*u0^2;

UU = ones([3, length(xx)]);

% density
UU(1,:) = [rho0*g_x(1), 2*gam*p_i/(gam-1).*g_x(2:end)];
% UU(1,:) = rho0.*g_x;

% velocityx
UU(2,:) = u0.*UU(1,:);
% UU(2,:) = u0.*UU(1,:);

%Energy
UU(3,:) = [E0*UU(1,1), E_i.*UU(1,2:end)];
% UU(3,:) = E0.*UU(1,:);

% UU(:,end,:) = UUexit.*g_x(end).*[0;1;1].*UUexit(1);

UU = repmat(UU,1,1,3);

%% Setup Time Simulation
time = 0;
tEnd = 10*5;

%% Run Simulation

while length(time) < 20
    
    % Update Field Values and boundary conditions
    UU(:,:,1:2) = UU(:,:,2:3);
    
    % Outflow Boundary Condition
    UU(:,end,2) = 4/3.*UU(:,end-1,2) - 1/3.*UU(:,end-2,2);
    
    % Update Flux
%     PP = (gam - 1).*(UU(:,3,2) - 0.5.*UU(:,2,2).^2./UU(:,1,2));
    FF =  flux(UU./g_x);
    
    CFL1 = dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    CFL2 = 2.*(UU(2,2:end-1,end)./UU(1,2:end-1,end)).*dt./(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    cflFactor = 1;
    if (~isempty(CFL1) && ((max(abs(CFL1(:))) > cfl)||(max(abs(CFL2(:))) ~=cfl))) && any(UU(2,:,2) > 0)
        cflFactor = min(cfl/max(abs(CFL2(:))),cfl/max(abs(CFL2(:))));
        dt = dt.*cflFactor;
        
        if abs(log10(cflFactor)) > 0.1
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps\n');
            fprintf('New time step:%0.5e\n', dt);
        end
        
    end    
    
    % Calculate Pressure
    PP = (gam-1).*(UU(3,:,2) - 0.5.*(UU(2,:,2).^2) / UU(1,:,2))./g_x;
%     PP(end) = p_i;
    
    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (UU(:,2:end-1,1).*(1 - CFL1) - dt.*(FF(:,3:end) - FF(:,1:end-2))./dx + [0;1;0].*PP(2:end-1).*dgdx(2:end-1) + 2.*dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx).*UU(:,3:end,2) + epsFunc(UU(:,1:end-1,2), dx).*UU(:,1:end-2,2)))./(1 + CFL1);
    time(end+1) = time(end)+dt;
    
    figure(1);
    for i = 1:size(UU,1)
        if i == 1
            plot(xx, UU(i,:,end)./g_x);
        else
            plot(xx, UU(i,:,end)./UU(1,:,end));
        end
        hold on;
    end

    PP = (gam-1).*(UU(3,:,2) - 0.5.*(UU(2,:,2).^2) / UU(1,:,2))./g_x;
    plot(xx,PP);

    hold off;
%     title('Sod Shock Tube');
    legend('\rho', 'u', 'e', 'p', 'Location', 'Best');
    drawnow;
    
end

%% Post Process

% UU = UU./g_x;

figure(1);
for i = 1:size(UU,1)
    if i == 1
        plot(xx, UU(i,:,end)./g_x);
    else
        plot(xx, UU(i,:,end)./UU(1,:,end));
    end
    hold on;
end

PP = (gam-1).*(UU(3,:,2) - 0.5.*(UU(2,:,2).^2) / UU(1,:,2))./g_x;
plot(xx,PP);

hold off;
title('Sod Shock Tube');
legend('\rho', 'u', 'e', 'p', 'Location', 'Best');