clc;
clear;
close all;

dirName = 'isentropicNozzle';

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
dt = 0.001;

%% Nozzle Geometry
g_x = 1 + (2.*xx-1).^2;
dgdx = [4*(2*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];

%% Fluid Properties
gam_list = 1.4:.1:2;
cfl = 1;

%% Perform Runs
for ii = 1:1%length(gam_list)
    % Assign Gamma
    gam = gam_list(ii);
    
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
        UU(:,end,:) = 5/2.*UU(:,end-1,:) - 2.*UU(:,end-2,:) + 0.5.*UU(:,end-3,:); % - 1/3.*UU(:,end-2,2:3);

        % 2) Fix Rho
        UU(1,end,:) = g_x(end).*((gam.*p_i).^(1/gam));
        
        % Check Stability?
        epsE = epsFunc(UU(:,2:end,2), dx);
        epsW = epsFunc(UU(:,1:end-1,2), dx);
        CFL = 0.5*dt - 0.5.*(epsE + epsW).*(dt/dx).^2;
        if visc_t >= max(CFL(:))
            visc_t = max(CFL(:))./1.5;
            fprintf('Time viscosity changed: %0.5f\n', visc_t);
        end            
        
        % Calculate Next Time-Step
        aa = reshape([[0;0], 0.5.*epsW./(dx^2), [0;0]]',[],1);
        bb = reshape([[1;1], 0.5.*(-(epsE + epsW)./dx^2) - visc_t./dt^2 - 1./(2*dt), [1;1]]',[],1);
        cc = reshape([[0;0], 0.5.*epsE./(dx^2), [0;0]]',[],1);
        dd = reshape([UU(:,1,3),    -2.*visc_t.*UU(:,2:end-1,2)./dt^2 ...
                                    + (FF(:,3:end) - FF(:,1:end-2))./(2*dx)...
                                    - [0;1].*PP(2:end-1).*dgdx(2:end-1)...
                                    + (visc_t./dt^2-0.5/dt+0.5.*(epsE + epsW)./dx^2).*UU(:,2:end-1,1)...
                                    - 0.5.*(epsE.*UU(:,3:end,1)+epsW.*UU(:,1:end-2,1))./dx^2,...
                                    UU(:,end,3)]',[],1);
        UU(:,:,3) = reshape(thomas3(aa,bb,cc,dd),size(UU,2),2)';
        time(end+1) = time(end)+dt;

        % Update Flux and Pressure
        FF =  fluxFuncIsentropic(UU(:,:,3)./g_x, gam).*g_x;
        PP = ((UU(1,:,3)./g_x).^gam)./gam;

        res(end+1,:) = reshape(max(abs((epsFunc(UU(:,2:end,3), dx).*(UU(:,3:end,3)-UU(:,2:end-1,3))-epsFunc(UU(:,1:end-1,3), dx).*(UU(:,2:end-1,3)-UU(:,1:end-2,3)))./(dx^2)...
            - ((FF(:,3:end) - FF(:,1:end-2))./(2.*dx) - [0;1].*PP(2:end-1).*dgdx(2:end-1))), [], 2), 1, 2);
    
        figure(1);
        semilogy(res);
        legend('\rho', 'u', 'Location', 'BestOutside');
%         for i = 1:size(UU,1)
%             if i == 1
%                 plot(xx, UU(i,:,end)./g_x, 'o-');
%             else
%                 plot(xx, UU(i,:,end)./UU(1,:,end), 'o-');
%             end
%             hold on;
%         end
% 
%         plot(xx,PP);
%         legend('\rho', 'u', 'p', 'Location', 'BestOutside');
        title(['Run: ' num2str(ii)]);
        drawnow;
        hold off;

    end
    
%% Post Process
    explicit = load([dirName '\gam_1']);
    figure(1);
    semilogy(res);
    legend('\rho', 'u', 'Location', 'BestOutside');
    title('Implicit Scheme');
    set(gca, 'FontSize', 14);
    saveas(gcf, 'implicit_residual', 'pdf');
    
    figure(2);
    semilogy(explicit.res);
    legend('\rho', 'u', 'Location', 'BestOutside');
    title('Explicit Scheme');
    set(gca, 'FontSize', 14);
    saveas(gcf, 'explicit_residual', 'pdf');

end