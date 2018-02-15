clc;
clear;
close all;

%% Generate Grid

% dx = 0.0025;
% xx = (0:dx:1)';
xx = linspace(0,1,151);
dx = xx(2) - xx(1);
dt = 0.1;

%% Nozzle Geometry
ind1 = xx>=0 & xx <=0.5;
ind2 = xx>0.5;
alpha = 0.4;
% g_x = [9/4.*(2.*xx(ind1) - 1).^4 - 1.5.*(2.*xx(ind1)-1).^2 + 5/4, (1-alpha).*(9/4.*(2.*xx(ind2)-1).^4 - 1.5.*(2.*xx(ind2)-1).^2) + 5/4];
% dgdx = [0, (g_x(3:end) - g_x(1:end-2))./(2.*dx), 0];

g_x = 1 + (2.*xx-1).^2;
dgdx = 4.*(2.*xx-1);

%% Fluid Properties
gam = 1.4;
cfl = 0.5;

flux =@(U) [ U(2,:,2);...
            0.5.*(3-gam).*(U(2,:,2).^2)./U(1,:,2) + (gam-1).*U(3,:,2);...
            -0.5.*(gam-1).*(U(2,:,2).^3)./(U(1,:,2).^2) + gam.*U(2,:,2).*U(3,:,2)./U(1,:,2)] .* g_x;

% flux =@(U) [ U(2,:,2);...
%             0.5.*(3-gam).*(U(2,:,2).^2)./U(1,:,2) + (gam-1).*E0;...
%             ].*g_x;

%% Initialize Field Values

% Inlet Conditions -> s0 = 0 (isentropic relations)
rho0 = 1.5;%(0.5.*(gam+1))^(1/(gam-1));
u0 = 1/3;
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;

% Exit Conditions
p_i = 1;

% UU = ones([2, length(xx)]);

% UU = [rho0; rho0*u0; rho0*E0].*g_x;
UU = [rho0.*g_x; rho0.*linspace(u0, 0.4, length(g_x)).*g_x; rho0.*E0.*g_x];
% UU(:,end,:) = [rho_i; rho_i.*u_i; rho_i.*E_i] .*g_x(end);
UU = repmat(UU,1,1,3);


%% Setup Time Simulation?
time = 0;
tEnd = 10*5;
res = ones(size(reshape(max(UU(:,:,3) - UU(:,:,2), [], 2)./dt, 1,3)));

%% Run Simulation

while length(time) < 2e3 || norm(res(end,:)) > 1e-3
    
    % Update Field Values and boundary conditions
    UU(:,:,1:2) = UU(:,:,2:3);
    
    % Outflow Boundary Condition
    % 1) Extrapolate rho and U
%     UU(:,end,2:3) = 4/3.*UU(:,end-1,2:3) - 1/3.*UU(:,end-2,2:3);
    UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);
    
    % 2) Calculate E from extrapolation
    UU(3,end,2:3) = (p_i.*g_x(end)/(gam-1) + 0.5.*(UU(2,end,2:3).^2)./UU(1,end,2:3));
    
    % Update Flux
    FF =  flux(UU./g_x);
    
    % Calculate CFL
    CFL1 = dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    CFL2 = (UU(2,2:end-1,end)./UU(1,2:end-1,end)).*dt./(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
    cflFactor = 1;
    if (~isempty(CFL1) && ((max(abs(CFL1(:))) > cfl)||(max(abs(CFL2(:))) ~=cfl))) && any(UU(2,:,2) > 0)
        cflFactor = min(cfl/max(abs(CFL2(:))),cfl/max(abs(CFL2(:))));
        dt = dt.*cflFactor;
        CFL1 = CFL1 .* (1 + cflFactor)/2;
        CFL2 = CFL2 .* (1+ cflFactor)/2;
        
        if abs(log10(cflFactor)) > 0.1
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps\n');
            fprintf('New time step:%0.5e\n', dt);
        end
        
    end    
    
    % Calculate Pressure
    PP = (gam-1).*(UU(3,:,2) - 0.5.*(UU(2,:,2).^2) ./ UU(1,:,2))./g_x;
    
    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (UU(:,2:end-1,1).*(1 - CFL1) + (1+1/cflFactor).*dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx).*UU(:,3:end,2) + epsFunc(UU(:,1:end-1,2), dx).*UU(:,1:end-2,2)) - (1+1/cflFactor).*dt.*((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1;0].*PP(2:end-1).*dgdx(2:end-1)))./(1 + CFL1);
    time(end+1) = time(end)+dt;
    
    res(end+1,:) = reshape(max(UU(:,:,3) - UU(:,:,2), [], 2)./dt, 1,3);
    
    figure(1);
%     semilogy(res);
    for i = 1:size(UU,1)
        if i == 1
            plot(xx, UU(i,:,end)./g_x);
        else
            plot(xx, UU(i,:,end)./UU(1,:,end));
        end
        hold on;
    end
    
    plot(xx,PP);

    hold off; 
%     title('Sod Shock Tube');
    legend('\rho', 'u', 'e', 'p', 'Location', 'BestOutside');
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
% PP = (gam-1).*(UU(1,:,2).*E0 - 0.5.*(UU(2,:,2).^2) ./ UU(1,:,2))./g_x;
plot(xx,PP);

hold off;
title('Sod Shock Tube');
legend('\rho', 'u', 'e', 'p', 'Location', 'Best');