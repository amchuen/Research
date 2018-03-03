clc;
clear;
close all;

dirName = 'testStability';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');
addpath('../../Matrix Solvers/');

%% Generate Grid

xx = linspace(0.5,1,151);
dx = xx(2) - xx(1);

visc_x = 0.002;
visc_t = 5e-5;
dt = 0.0025;
ratio_list = [1.1, 1.01, 1.001, 0.999];

%% Nozzle Geometry
g_x = 1 + (2.*xx-1).^2;
dgdx = [4*(2*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];

%% Fluid Properties
% gam_list = 1.4:.1:2;
gam = 2;
cfl = 1;

%% Perform Runs
for ii = 1:1%length(ratio_list)
    % Assign Gamma
%     gam = gam_list(ii);
    ratio = ratio_list(ii);
    
    % Initialize Field Values
    % Inlet Conditions -> s0 = 0 (isentropic relations)
    if xx(1) == 0
        rho0 = 1.5;%(0.5.*(gam+1))^(1/(gam-1));
        u0 = 1/3;
    elseif xx(1) == 0.5
        rho0 = 1;
        u0 = 1;
    end
    p0 = (rho0^gam)/gam;
    E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;
    H0 = 3;

    % Exit Conditions
    p_i = 1;
    % rho_e = 1.20135; %1.2013512490310365;
    u_e = 0.416198;% 0.4162064994502009;
    rho_e = 1/(g_x(end)*u_e);

    % UU = [rho0; rho0*u0; rho0*E0].*g_x;
    UU = [rho0.*g_x; rho0.*linspace(u0, 0.4, length(g_x)).*g_x];
    UU(:,end) = [rho_e; rho_e*u_e].*g_x(end);
    UU = repmat(UU,1,1,3);

    
    % Setup Time Simulation?
    time = 0;
    tEnd = 10*5;
    FF =  fluxFuncIsentropic(UU(:,:,3)./g_x, gam).*g_x;
    PP = ((UU(1,:,3)./g_x).^gam)./gam;
    res = reshape(max(abs((visc_x.*(UU(:,3:end,3)-UU(:,2:end-1,3))-visc_x.*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
            - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1].*PP(2:end-1).*dgdx(2:end-1))), [], 2), 1, 2);
    dtLast = dt;
        
    while length(time) < 3 || norm(res(end,:)) > 1e-3

       % Update Field Values and boundary conditions
        UU(:,:,1:2) = UU(:,:,2:3);

        % Outflow Boundary Condition
        % 1) Extrapolate rho and E
        UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);

        % 2) Fix Rho
        UU(1,end,2:3) = g_x(end).*((gam.*p_i).^(1/gam));

        % Calculate CFL
        visc_x = (epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
        CFL1 = dt./(dx^2).*(epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx));
        U_0 = abs(UU(2,2:end-1,end)./UU(1,2:end-1,end));
        U_pA = abs(U_0 + sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
        U_mA = abs(U_0 - sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
        Umax = max([max(U_0(:)), max(U_pA(:)), max(U_mA(:))]);
        dt = dx./(sqrt(max(visc_x(:))/visc_t - Umax^2));
        if visc_t <= dt.*max(CFL1(:))./2
            visc_t = ratio.*dt.*max(CFL1(:))./2;
            fprintf('Changing time viscosity\nNew time visc: %0.5f\n', visc_t);
        end
        
        % Calculate Next Time-Step
        UU(:,2:end-1,3) = ((epsFunc(UU(:,2:end,2), dx).*UU(:,3:end,2) + epsFunc(UU(:,1:end-1,2), dx).*UU(:,1:end-2,2))./(dx^2)...
                            - UU(:,2:end-1,2).*((epsFunc(UU(:,2:end,2), dx) + epsFunc(UU(:,1:end-1,2), dx))./(dx^2) - 2.*visc_t./(dt^2))...
                            - (FF(:,3:end) - FF(:,1:end-2))./(2*dx) + [0;1].*PP(2:end-1).*dgdx(2:end-1)...
                            - UU(:,2:end-1,1).*(visc_t./(dt^2)-0.5/dt))./(visc_t./(dt^2) + 0.5/dt);
        time(end+1) = time(end)+dt;

        % Update Flux and Pressure
        FF =  fluxFuncIsentropic(UU(:,:,3)./g_x, gam).*g_x;
        PP = ((UU(1,:,3)./g_x).^gam)./gam;

        res(end+1,:) = reshape(max(abs((epsFunc(UU(:,2:end,3), dx).*(UU(:,3:end,3)-UU(:,2:end-1,3))-epsFunc(UU(:,1:end-1,3), dx).*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
            - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1].*PP(2:end-1).*dgdx(2:end-1))), [], 2), 1, 2);
    
        figure(1);
        semilogy(res);
        title(['Run: ' num2str(ii)]);
        legend('\rho', 'u', 'Location', 'BestOutside');
        drawnow;

    end
    
%% Post Process
    fileName = ['ratio_' num2str(ii)];
    save([dirName '\' fileName]);
    
    figure();
    plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
    plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
    plot(xx, PP, '^-');
    title(['Ratio = ' num2str(ratio)]);
    legend('\rho', 'u', 'P', 'Location', 'Best');
    drawnow;
    saveas(gcf, [dirName '\' fileName]);
    
    close all;

end