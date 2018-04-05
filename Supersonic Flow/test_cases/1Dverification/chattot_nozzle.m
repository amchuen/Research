clc;
clear;
close all;

dirName = 'chattot_nozzle';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');
addpath('../../Matrix Solvers/');

%% Generate Grid

xx = linspace(0,1,151);
dx = xx(2) - xx(1);

visc_x = 0.005;
% beta = 5e-6;
ratio_list = linspace(0.975, 1, 7);% [0.85, 0.9, 0.999, 1];
tol = 1e-3;

%% Nozzle Geometry
aa = 0.4;
g_x = [ 9/4.*(2*xx(xx<0.5)-1).^4 - 1.5.*(2*xx(xx<0.5)-1).^2 + 5/4,...
        (1-aa).*(9/4.*(2*xx(xx>=0.5)-1).^4 - 1.5.*(2*xx(xx>=0.5)-1).^2)+5/4];
dgdx = [2.*9.*(2.*xx(1)-1).^3 - 3.*2.*(2.*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), (1-aa).*(2.*9.*(2.*xx(end)-1).^3 - 3.*2.*(2.*xx(end)-1))];

figure();plot(xx, g_x);

%% Fluid Properties
% gam_list = 1.4:.1:2;
gam = 2;
cfl = 1;

%% Perform Runs
for ii = 8:(7+length(ratio_list))
    % Assign Gamma
%     gam = gam_list(ii);
    ratio = ratio_list(ii-7);
    dt = dx^2;
    beta = visc_x*(dt^2)/(ratio*dx^2);
    
    % Initialize Field Values
    % Inlet Conditions -> s0 = 0 (isentropic relations)
    if xx(1) == 0
        rho0 = (0.5.*(gam+1))^(1/(gam-1));
        u0 = 0;
    elseif xx(1) == 0.5
        rho0 = 1;
        u0 = 1;
    end
    p0 = (rho0^gam)/gam;
    E0 = p0/((gam-1)*rho0) + 0.5.*u0^2;
    H0 = 3;

    % Exit Conditions
    p_i = 0.9;
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
    UU_x = UU(:,xx==0.65,3);
    dtLast = dt;
    
    figure(1);
    resRho = semilogy(res(:,1)); hold on;
    resU = semilogy(res(:,2));
    legend('\rho', 'u', 'Location', 'BestOutside');
    title('Residual');
    movegui(gcf, 'west');
    
%     figure(2);
%     Rho_x = plot(UU_x(1,:)); hold on;
%     U_x = plot(UU_x(2,:)./UU_x(1,:));
%     legend('\rho', 'u', 'Location', 'BestOutside');
%     title('Oscillation');
%     movegui(gcf, 'east');
%     Ux_fig = figure();
        
    while length(time) < 3 || norm(res(end,:)) > tol

       % Update Field Values and boundary conditions
        UU(:,:,1:2) = UU(:,:,2:3);

        % Outflow Boundary Condition
        % 1) Extrapolate rho and E
        UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);

        % 2) Fix Rho
        UU(1,end,2:3) = g_x(end).*((2.*gam.*p_i).^(1/gam));
        
        % Check CFL
        U_0 = abs(UU(2,2:end-1,end)./UU(1,2:end-1,end));
        U_pA = abs(U_0 + sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
        U_mA = abs(U_0 - sqrt(gam.*PP(2:end-1).*g_x(2:end-1)./UU(1,2:end-1,end)));%.*g_x(2:end-1));
        Umax = max([max(U_0(:)), max(U_pA(:)), max(U_mA(:))]);
        if Umax*dt/dx ~= 1 %Umax*dt/dx > 1
%             error('CFL condition not met!\n');
            dt = dx./Umax;
            beta = visc_x*(dt^2)/(ratio*dx^2);
            if abs(log10(dt/dtLast)) > 0.301
                dtLast = dt;
                fprintf('Time-step changing!\nNew time step: %0.5e\nNew beta:%0.5e\n\n', dt, beta);
            end
        end

        % Calculate Next Time-Step
        UU(:,2:end-1,3) = ((UU(:,3:end,2) + UU(:,1:end-2,2)).*(visc_x.*dt^2)./(beta.*dx^2) + 2.*UU(:,2:end-1,2).*(1-visc_x*dt^2/(beta*dx^2))-(1-dt/(2*beta))*UU(:,2:end-1,1)...
                            - (dt^2./beta).*((FF(:,3:end) - FF(:,1:end-2))./(2*dx) - [0;1].*PP(2:end-1).*dgdx(2:end-1)))./(1+0.5*dt/beta);
        time(end+1) = time(end)+dt;

        % Update Flux and Pressure
        FF =  fluxFuncIsentropic(UU(:,:,3)./g_x, gam).*g_x;
        PP = ((UU(1,:,3)./g_x).^gam)./gam;

        res(end+1,:) = reshape(max(abs((visc_x.*(UU(:,3:end,3)-UU(:,2:end-1,3))-visc_x.*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
            - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1].*PP(2:end-1).*dgdx(2:end-1))), [], 2), 1, 2);
%         UU_x(:,end+1) = UU(:,xx==0.65,3);
        
        % Plot Residuals
        resRho.YData(end+1) = res(end,1);
        resU.YData(end+1) = res(end,2);
%         Rho_x.YData(end+1) = UU_x(1,end);
%         U_x.YData(end+1) = UU_x(2,end)./UU_x(1,end);
        drawnow;

    end
    
%% Post Process
    fileName = ['ratio_' num2str(ii)];
    save([dirName '\' fileName]);
    
    figure(1);
    saveas(gcf, [dirName '\resFig_' fileName]);
    
    figure(2);
    saveas(gcf, [dirName '\oscFig_' fileName]);
    
    figure();
    plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
    plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
    plot(xx, PP, '^-');
    title(['Ratio = ' num2str(ratio)]);
    legend('\rho', 'u', 'P', 'Location', 'Best');
    drawnow;
    saveas(gcf, [dirName '\' fileName]);
    
    fprintf('Simulation complete!\nNo. of iterations: %i\n', length(res));
    if norm(res(end)) < tol
        fprintf('Success: Yes\n\n');
    else
        fprintf('Success: No\n\n');
    end
    
    close all;

end