clc;
clear;
close all;

%% Generate Grid

% dx = 0.0025;
% xx = (0:dx:1)';
xx = linspace(0.5,1,101);
dx = xx(2) - xx(1);
dt = 5*dx^2;

%% Nozzle Geometry
g_x = 1 + (2.*xx-1).^2;
% Extend nozzle
xx(end+1:end+2) = [xx(end)+dx, xx(end)+2*dx];
g_x(end+1:end+2) = g_x(end);

dgdx = [4*(2*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];

%% Fluid Properties
gam = 1.4;
cfl = 0.125;

%% Initialize Field Values

% Inlet Conditions -> s0 = 0 (isentropic relations)
% rho0 = 1.5;%(0.5.*(gam+1))^(1/(gam-1));
% u0 = 1/3;
rho0 = 1;
u0 = 1;
p0 = (rho0^gam)/gam;
E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;

% Exit Conditions
p_i = 1;
rho_e = 1.2013512490310365;
u_e = 0.4162064994502009;


% UU = ones([2, length(xx)]);

% UU = [rho0; rho0*u0; rho0*E0].*g_x;
UU = [rho0.*g_x; rho0.*linspace(u0, 0.4, length(g_x)).*g_x; rho0.*E0.*g_x];
% UU(:,end,:) = [rho_i; rho_i.*u_i; rho_i.*E_i] .*g_x(end);
UU(:,end) = [rho_e; rho_e*u_e; rho_e*E0].*g_x(end);
UU = repmat(UU,1,1,3);


%% Setup Time Simulation?
time = 0;
tEnd = 10*5;
FF =  fluxFunc(UU(:,:,3), gam);
PP = (gam-1).*(UU(3,:,3) - 0.5.*(UU(2,:,3).^2) ./ UU(1,:,3))./g_x;
res = reshape(max((epsFunc(UU(:,2:end,3), dx).*(UU(:,3:end,3)-UU(:,2:end-1,3))-epsFunc(UU(:,1:end-1,3), dx).*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
        - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1;0].*PP(2:end-1).*dgdx(2:end-1)), [], 2), 1, 3);
dtLast = dt;

maxRes = max(res, [], 1);

%% Run Simulation

while norm(res(end,:)./maxRes) > 1e-3
    
    % Update Field Values and boundary conditions
    UU(:,:,1:2) = UU(:,:,2:3);
    
    % Outflow Boundary Condition
    % 1) Extrapolate rho and E
    UU(:,end,2:3) = 4/3.*UU(:,end-1,2:3) - 1/3.*UU(:,end-2,2:3);
%     UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);
    
    % 2) Fix U
%     UU(2,end,2:3) = UU(1,end,2:3).*u_e;
    UU(3,end,2:3) = (p_i.*g_x(end)/(gam-1) + 0.5.*(UU(2,end,2:3).^2)./UU(1,end,2:3));
    
    % Calculate CFL
    CFL1 = 0.5.*dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
%     U_0 = abs(UU(2,2:end-1,end)./UU(1,2:end-1,end));
%     U_pA = abs(U_0 + sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
%     U_mA = abs(U_0 - sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
%     Umax = max([max(U_0(:)), max(U_pA(:)), max(U_mA(:))]);
%     CFL2 = Umax.*dt./(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
%     CFL3 = (UU(2,2:end-1,end)./UU(1,2:end-1,end)).*dt./dx;
    cflFactor = 1;
    if (max(abs(CFL1(:))) ~= cfl) %(~isempty(CFL1) && ((max(abs(CFL1(:))) > cfl)||(max(abs(CFL2(:))) ~=cfl))) && any(UU(2,:,2) > 0)
%         cflFactor = min([cfl/max(abs(CFL2(:))),cfl/max(abs(CFL1(:)))]);
        cflFactor = cfl/max(abs(CFL1(:)));
        dt = dt.*cflFactor;
        CFL1 = CFL1 .* cflFactor;
%         CFL2 = CFL2 .* cflFactor;
        
        if abs(log10(dtLast/dt)) > 0.1
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps\n');
            fprintf('New time step:%0.5e\n', dt);
            dtLast = dt;
        end
        
    end    
    
    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (UU(:,2:end-1,2).*(1 - 2.*CFL1)...
                        + dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx).*UU(:,3:end,2) + epsFunc(UU(:,1:end-1,2), dx).*UU(:,1:end-2,2))...
                        - dt.*((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1;0].*PP(2:end-1).*dgdx(2:end-1)));
    time(end+1) = time(end)+dt;
    
    % Update Flux and Pressure
    FF =  fluxFunc(UU(:,:,3), gam);
    PP = (gam-1).*(UU(3,:,3) - 0.5.*(UU(2,:,3).^2) ./ UU(1,:,3))./g_x;
    
    res(end+1,:) = reshape(max(abs((epsFunc(UU(:,2:end,3), dx).*(UU(:,3:end,3)-UU(:,2:end-1,3))-epsFunc(UU(:,1:end-1,3), dx).*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
        - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1;0].*PP(2:end-1).*dgdx(2:end-1))), [], 2), 1, 3);
    if size(res,1) <= 10
        maxRes = max(res, [], 1);
    else
        maxRes = max(res(1:10,:), [], 1);
    end
    
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