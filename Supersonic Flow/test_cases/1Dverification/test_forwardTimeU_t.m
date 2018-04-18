clc;
clear;
close all;

dirName = 'testForwardTimeU_t';

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

addpath(dirName);
addpath('fluxSchemes\');
addpath('viscositySchemes\');
addpath('../../../Matrix Solvers/');

% Notes:
    % - Test convergence of different time discretization using different
    %   grid size and different time-diffusion ratios
    % - Results are sorted by time-discretization method, then by ratio,
    %   then finally by grid size

%% Independent Variables
% vary ratio, grid density, and time-discretization
ratio_list = [linspace(0.85, 0.95, 7), linspace(0.975,1, 7)];% [0.85, 0.9, 0.999, 1];
nPts_list = [101, 201, 401];
tol = 1e-3;
maxIter = 0.5e5;

%% Fluid Properties and Boundary Conditions
gam = 2;
cfl = 1;

% Inlet Conditions -> s0 = 0 (isentropic relations)
rho0 = 1;
u0 = 1;
p0 = (rho0^gam)/gam;
H0 = 3;

% Exit Conditions -> give one, extrapolate the rest
p_i = 1;

for i = 1:2
    %% Set Time Discretization
    if i == 1
        folder= 'forward_time';
    elseif i == 2
        folder = 'central_time';
    end
    
    % Collect data for n residuals vs ratio for three different grid sizes
    x_n1 = [];
    y_n1 = [];
    x_n2 = [];
    y_n2 = [];
    x_n3 = [];
    y_n3 = [];
    
    for ii = 1:length(ratio_list)
        %% Set Ratio
        ratio = ratio_list(ii);
        subFolder = ['ratio_' num2str(ii)];
        
        for iii = 1:length(nPts_list)
            %% Set Grid Discretization
            xx = linspace(0.5,1,nPts_list(iii));
            dx = xx(2) - xx(1);
            dt = 0.4*dx;

            g_x = 1 + (2.*xx-1).^2;
            dgdx = [4*(2*xx(1)-1), (g_x(3:end)-g_x(1:end-2))./(2*dx), 0];
            
            subSubFolder = num2str(nPts_list(iii));

            %% Perform Runs
            % Initialize Field Values
            UU = [rho0.*g_x; rho0.*linspace(u0, 0.4, length(g_x)).*g_x];
%             UU(:,end) = [rho_e; rho_e*u_e].*g_x(end);
            UU = repmat(UU,1,1,3);

            % Setup Time Simulation
            time = 0;
            tEnd = 10*5;
            [flux, Umax, ~, visc_beta] = fx_2Diff(@fluxFuncIsentropic, @epsConst, UU(:,:,end), g_x, gam, dx);
            res = reshape(max(abs(flux), [], 2), 1, 2);
            UU_x = UU(:,xx==0.65,3);
            dtLast = dt;
            beta = visc_beta.*(dt^2)/(ratio*dx^2);

            figure(1);
            resRho = semilogy(res(:,1)); hold on;
            resU = semilogy(res(:,2));
            legend('\rho', 'u', 'Location', 'BestOutside');
            title('Residual');
            movegui(gcf, 'west');

            figure(2);
            Rho_x = plot(UU_x(1,:)); hold on;
            U_x = plot(UU_x(2,:)./UU_x(1,:));
            legend('\rho', 'u', 'Location', 'BestOutside');
            title('Oscillation');
            movegui(gcf, 'east');
            %     Ux_fig = figure();

            while length(time) < maxIter && norm(res(end,:)) > tol

               % Step-forward Field Values
                UU(:,:,1:2) = UU(:,:,2:3);

                % Check CFL
                if Umax*dt/dx > cfl
                    dt = cfl.*dx./Umax;
                    if abs(log10(dt/dtLast)) > 1e-4
                        dtLast = dt;
                        fprintf('Time-step changing!\nNew time step: %0.5e\n', dt);
                    end
                end

                % Calculate Next Time-Step
                if i == 1
                    beta = 0.5.*(2.*visc_beta.*(dt^2)/(ratio*dx^2) - dt);
                    UU(:,2:end-1,3) = ((2.*beta./(dt^2) + 1/dt).*UU(:,2:end-1,2)-(beta./(dt^2)).*UU(:,2:end-1,1)-flux)./(beta./(dt^2)+1/dt);
                elseif i == 2
                    beta = visc_beta.*(dt^2)/(ratio*dx^2);
                    UU(:,2:end-1,3) = ((2.*beta./(dt^2)).*UU(:,2:end-1,2)-(beta./(dt^2)-0.5./dt).*UU(:,2:end-1,1)-flux)./(beta./(dt^2)+0.5/dt);
                end
                time(end+1) = time(end)+dt;

                % Update Outflow Boundary Condition
                % 1) Extrapolate rho and E
                UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);

                % 2) Fix Rho
                UU(1,end,2:3) = g_x(end).*((gam.*p_i).^(1/gam));

                % Update Flux and Pressure
                [flux, Umax, ~, visc_beta] = fx_2Diff(@fluxFuncIsentropic, @epsConst, UU(:,:,end), g_x, gam, dx);
                res(end+1,:) = reshape(max(abs(flux), [], 2), 1, 2);
                UU_x(:,end+1) = UU(:,xx==0.65,3);

                % Plot Residuals
                resRho.YData(end+1) = res(end,1);
                resU.YData(end+1) = res(end,2);
                Rho_x.YData(end+1) = UU_x(1,end);
                U_x.YData(end+1) = UU_x(2,end)./UU_x(1,end);
%                 drawnow;

            end

            %% Post Process
            fprintf('Simulation complete!\nNo. of iterations: %i\n', length(res));
            fprintf('Time discretization: %s\n', folder);
            if norm(res(end)) < tol
                fprintf('Success: Yes\n\n');
                folderName = [dirName '\' folder '\' subFolder '\' subSubFolder];
                if ~exist(folderName, 'dir')
                    mkdir(folderName);
                end

                save([folderName '\output']);

                figure(1);
                saveas(gcf, [folderName '\residual']);

                figure(2);
                saveas(gcf, [folderName '\oscillation']);

                figure();
                [~, PP] = fluxFuncIsentropic(UU(:,:,end)./g_x, gam);
                plot(xx, UU(1,:,3)./g_x, '*-'); hold on;
                plot(xx, UU(2,:,3)./UU(1,:,3), 'o-');
                plot(xx, PP, '^-');
                title(['Ratio = ' num2str(ratio)]);
                legend('\rho', 'u', 'P', 'Location', 'Best');
                drawnow;
                saveas(gcf, [folderName '\profile']);
                
                % Add to Data plot
                if iii == 1
                    x_n1(end+1) = ratio;
                    y_n1(end+1) = length(res);
                elseif iii == 2
                    x_n2(end+1) = ratio;
                    y_n2(end+1) = length(res);
                elseif iii == 3
                    x_n3(end+1) = ratio;
                    y_n3(end+1) = length(res);
                end
                
            else
                fprintf('Success: No\n\n');
            end
            
            close all;
        end
    end
    
    % Plot data
    figure();
    semilogy(x_n1, y_n1, 'o-');
    hold on;
    semilogy(x_n2, y_n2, '*-');
    semilogy(x_n3, y_n3, '^-');
    xlabel('Ratio of viscosities');
    ylabel('No. of iterations.');
    if i == 1
        title('Forward Time U_t');
    elseif i == 2
        title('Central Time U_t');
    end
    legend('101 pts', '201 pts', '401 pts', 'Location', 'Best');
    saveas(gcf, [dirName '\' folder '\total_data']);
    
    close all;
end