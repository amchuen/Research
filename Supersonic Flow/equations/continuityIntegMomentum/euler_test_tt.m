clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dx = 0.025;
dy = 0.04;

% Field Axis Values
y_max = 10;
x_max = 10 + 20*dx;
x_min = -10-39*dx; %(-19*dx);
x_vals = x_min:dx:x_max;
y_vals = 0:dy:y_max;
[XX, YY] = meshgrid(x_vals, y_vals);

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.95;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_1 = 0.05;
eps_2 = eps_1;
tol = 1e-7;
min_tol = tol;
dt = dx^3;
three_level = 1; % switch between two-level and three-level scheme
CFL_on = 1;
iter_min = 800;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.1;
m_x = tand(8); % dy/dx
% x_vals = x_vals;
% dx = dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i+1) = (YY_B(i) - YY_B(i-1))/(2*dx);
end

BC.Vy_II = zeros(size(YY(end,:)));
BC.Vx_II = ones(size(XX(end,:)));
BC.Vx_I = ones(size(XX(:,1)));
BC.PHI_II = XX(end,:);
BC.PHI_I = XX(:,1);
BC.dyBdx = dyBdx;

%% START SOLUTION

PHI = zeros([size(XX), 3]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
PHI(:,:,3) = PHI(:,:,2);
RHO = ones([size(XX), 3]);
res = 1;

while (res(end) > tol)|| (length(res) < iter_min) % iterate through time
        
    % Reorganize three-level scheme
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI(:,:,3);
    
    RHO(:,:,1) = RHO(:,:,2);
    RHO(:,:,2) = RHO(:,:,3);
    
    if M0 ~= 0
        % Density Averages
        rho_e = [0.5.*(RHO(:,2:end,2) + RHO(:,1:(end-1),2)), 0.5.*(ones(size(RHO(:,end,2))) + RHO(:,end,2))];
        rho_w = [0.5.*(ones(size(RHO(:,1,2))) + RHO(:,1,2)), 0.5.*(RHO(:,1:(end-1),2) + RHO(:,2:end,2))];
        rho_n = [0.5.*(RHO(2:end,:,2) + RHO(1:(end-1),:,2)); 0.5.*(ones(size(RHO(end,:,2))) + RHO(end,:,2))];
        rho_s = 0.5.*([2.*RHO(1,:,2); RHO(2:end,:,2) + RHO(1:(end-1),:,2)]);
        
        % Density Gradient -> first order 
        rhoX_w = [(RHO(:,1,2) - ones(size(RHO(:,1,2))))./dx, (RHO(:, 2:end, 2) - RHO(:,1:(end-1),2))./dx];
        rhoX_e = [(RHO(:, 2:end, 2) - RHO(:,1:(end-1),2))./dx, (ones(size(RHO(:,end,2))) - RHO(:,end,2))./dx];
        rhoY_n = [(diff(RHO(:,:,2),1,1)./dy); (ones(size(RHO(end,:,2))) - RHO(end,:,2))./dy];
        rhoY_s = [zeros(size(RHO(1,:,2))); (diff(RHO(:,:,2),1,1)./dy)]; % enforce irrotationality at the body
        
        % Velocity -> first order
        phiX_w = [ones(size(PHI(:,1,2))), (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dx];
        phiX_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dx, ones(size(PHI(:,end,2)))];
        phiY_n = [(diff(PHI(:,:,2),1,1)./dy); (XX(end,:) - PHI(end,:,2))./dy];
        phiY_s = [dyBdx; (diff(PHI(:,:,2),1,1)./dy)]; % enforce irrotationality at the body

        phiX = 0.5.*(phiX_w + phiX_e);
        phiY = 0.5.*(phiY_n + phiY_s);
        q2_ij = phiX.^2 + phiY.^2;

        % Mass Flux, Steady Term
        mflux_ss = (rho_e.*phiX_e - rho_w.*phiX_w)./dx + (rho_n.*phiY_n - rho_s.*phiY_s)./dy;

        % Laplace of PHI
        phi_xx = (phiX_e - phiX_w)./dx;
        phi_yy = (phiY_n - phiY_s)./dy;
        phi_laplace = phi_xx + phi_yy;
        
         % Laplace of RHO
        rho_xx = (rhoX_e - rhoX_w)./dx;
        rho_yy = (rhoY_n - rhoY_s)./dy;
        rho_laplace = rho_xx + rho_yy;

        % Check CFL conditions
        CFL_i = (max(abs(phiX(:)))./dx + max(abs(phiY(:)))./dy)*dt;

        if CFL_i >= 1.0 && CFL_on
           fprintf('CFL condition not met!\n');
%            if CFL_on
           fprintf('Decreasing time steps!\n');
           dt = dt*0.8 / CFL_i;
           fprintf('New time step:%0.5f\n', dt);
%            end
%         elseif CFL_i < 0.8
%             dt = dt / CFL_i;
        end

        if three_level % three level scheme activated
            % Total Continuity Equation with unsteady term
            RHO(:,:,3) = ((2*eps_1./(dt^2)).*RHO(:,:,2) + (1/(2*dt) - eps_1/(dt^2)).*RHO(:,:,1) + eps_1.*rho_laplace - mflux_ss)./((eps_1/(dt^2) + 1/(2*dt)));
%             if any(abs((RHO(end,:,3)-1)) > 1e-2)
% %                 fprintf('Boundary condition not enforced for density!\n');
%                 RHO(end,abs((RHO(end,:,3)-1)) > 1e-5,3) = 1;
%             end
            
            % Total energy equation with unsteady term
            PHI(:,:,3) = ((2*eps_2./(dt^2)).*PHI(:,:,2) + (1/(2*dt) - eps_2/(dt^2)).*PHI(:,:,1) + eps_2.*phi_laplace - (RHO(:,:,2).^(gam-1) - 1)./((gam-1)*M0^2) - 0.5.*(q2_ij-1))./((eps_2/(dt^2) + 1/(2*dt)));
%             if any(abs((PHI(end,:,3)-XX(end,:)))./XX(end,:) > 1e-2)
% %                 fprintf('Boundary condition not enforced for potential!\n');
%                 PHI(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5,3) = XX(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5);
%             end
            
        else % use two-level time scheme
            RHO(:,:,3) = ((1/dt).*RHO(:,:,2) + eps_1.*rho_laplace - mflux_ss)./((1/(dt)));
            if any(abs((RHO(end,:,3)-1)) > 1e-5)
                fprintf('Boundary condition not enforced for density!\n');
                RHO(end,abs((RHO(end,:,3)-1)) > 1e-5,3) = 1;
            end

            % Total energy equation with unsteady term
            PHI(:,:,3) = ((1/dt).*PHI(:,:,2) + eps_2.*phi_laplace - (RHO(:,:,2).^(gam-1) - 1)./((gam-1)*M0^2) - 0.5.*(q2_ij-1))./((1/(dt)));
            if any(abs((PHI(end,:,3)-XX(end,:))) > 1e-5)
                fprintf('Boundary condition not enforced for potential!\n');
                PHI(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5,3) = XX(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5);
            end

        end

       % Checks
        if any(RHO(:,:,3) < 0)
            [r_ind, c_ind] = find(RHO(:,:,3) <0);
           fprintf('Check density!\n');
           fprintf('Issue at %0.2f, %0.2f\n', XX(r_ind,c_ind), YY(r_ind,c_ind));
           fprintf('Setting Negative Density to 0!\n');
           RHO(r_ind, c_ind,3) = 0;
        end

        if (any(abs(phiX(:)) > 2)) || (any(abs(phiY(:)) > 2))
            [r_X, c_X] = find(phiX > 2);
            [r_Y, c_Y] = find(phiY > 2);
           fprintf('Velocity too high!\n');
           fprintf('Issue at %0.2f, %0.2f\n', XX([r_X;r_Y],[c_X;c_Y]), YY([r_X;r_Y],[c_X; c_Y]));
           fprintf('Max Velocity: %0.5f\n\n', max(max(abs(phiX)), max(abs(phiY))));
        end

        if ~isreal(PHI(:,:,3))
           fprintf('Check phi!\n'); 
        end
    end
    
    % Run Residual Check
    difference = abs(PHI(:,:,3) - PHI(:,:,2));
    [res(end+1), ~] = max(difference(:));
    
    if (res(end) > 1e1) && (res(end-1) > 1e1)
        figure(1);semilogy(1:length(res), res);
        title('Residuals');
        xlabel('# of iterations');
        ylabel('Residual (Error)');
        error('Residual growing uncontrollably! Current residual at %0.5f', res(end)); 
    end

    % Initialize next time step
    if (length(res) > 50) && (mod(length(res), 50) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
    end
    
    % Break if residual starts to grow
    if ((res(end) > res(end-1)) && (res(end) < (min_tol))) && (length(res) > iter_min)
        fprintf('Iteration loop breaking at %i iterations.\n', length(res)) 
        break;
    end
    
    drawnow;
end

fprintf('\nCalculation Complete!\n');
%% POST PROCESS

% Calculate Velocities
% phiX = [ones(size(PHI(:,1,3))), (PHI(:,3:end,3) - PHI(:,1:(end-2),3))./(2*dx), ones(size(PHI(:,1,3)))];
% phiY = [dyBdx; (PHI(3:end,:,3) - PHI(1:(end-2),:,3))./(2*dy); zeros(size(dyBdx))];
% q2_ij = phiX.^2 + phiY.^2;

folderName = ['M_' num2str(M0)];

if ~exist([pwd '\euler_tt\' folderName], 'dir')
    mkdir([pwd '\euler_tt\' folderName]);
end

close all;
figure(1);semilogy(1:length(res), res);
title('Residual Plot');
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd '\euler_tt\' folderName '\residual_plot.png']);
saveas(gcf, [pwd '\euler_tt\' folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,round(RHO(:,:,3),5), 50)
title('Density (Normalized)');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_tt\' folderName '\density.png']);
saveas(gcf, [pwd '\euler_tt\' folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-(q2_ij), 50); %./((RR.*cos(TT)).^2)
title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_tt\' folderName '\cp_contour.png']);
saveas(gcf, [pwd '\euler_tt\' folderName '\cp_contour']);

% figure();
% % plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
% plot(fliplr(TT(1,:)).*180 ./ pi, 1-q2_ij(1,:));
% xlabel('\X');
% %     ylabel('\phi_{\theta}');
% ylabel('C_p');
% title('C_p on surface of Cylinder');
% % set(gca, 'Xdir', 'reverse');
% saveas(gcf, [pwd '\euler_tt\' folderName '\cp_surf.png']);
% saveas(gcf, [pwd '\euler_tt\' folderName '\cp_surf']);

figure(); % field potential
contourf(XX, YY, round(PHI(:,:,end),5), 50);
title('Field Potential, \Phi');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_tt\' folderName '\phi_pot.png']);
saveas(gcf, [pwd '\euler_tt\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, round(phiX,5), 50); %./((RR.*cos(TT)).^2)
title('\Phi_\theta velocity');
colorbar('eastoutside');
saveas(gcf, [pwd '\euler_tt\' folderName '\phi_theta.png']);
saveas(gcf, [pwd '\euler_tt\' folderName '\phi_theta']);

% figure(); % theta-dir velocity plots
% contourf(XX, YY, phiY, 50); %./((RR.*cos(TT)).^2)
% title('\Phi_R velocity');
% colorbar('eastoutside');
% saveas(gcf, [pwd '\euler_tt\' folderName '\phi_radius.png']);
% saveas(gcf, [pwd '\euler_tt\' folderName '\phi_radius']);

% figure();
% contourf(XX, YY, sqrt(M2_ij), 50);
% title('Mach Number');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\euler_tt\' folderName '\mach.png']);
% saveas(gcf, [pwd '\euler_tt\' folderName '\mach']);

%% Save Results
save([pwd '\euler_tt\' folderName '\results.mat']);

%% OLD CODE

%         for i = 2:(size(PHI,2)-1)
%             
%             if three_level % three level scheme activated
%                 % Total Continuity Equation with unsteady term
%                 RHO(1:(end-1),i,3) = ((2*eps_1./(dt^2)).*RHO(1:(end-1),i,2) + (1/(2*dt) - eps_1/(dt^2)).*RHO(1:(end-1),i,1) + eps_1.*rho_laplace - mass_con)./((eps_1/(dt^2) + 1/(2*dt)));
%                 RHO(end,i,3) = ones(size(RHO(end,i,3)));
% 
%                 % Total energy equation with unsteady term
%                 PHI(1:(end-1),i,3) = ((2*eps_2./(dt^2)).*PHI(1:(end-1),i,2) + (1/(2*dt) - eps_2/(dt^2)).*PHI(1:(end-1),i,1) + eps_2.*phi_laplace - (RHO(1:(end-1),i,2).^(gam-1) - 1)./((gam-1)*M0^2) - 0.5.*(q2_ij-1))./((eps_2/(dt^2) + 1/(2*dt)));
%                 PHI(end,i,3) = XX(end,i);
%             else % use two-level time scheme
%                 RHO(:,i,3) = ((1/dt).*RHO(:,i,2) + eps_1.*rho_laplace(:,i) - mass_con(:,i))./((1/(dt)));
%                 if abs((RHO(end,i,3)-1)) > 1e-5
%                     fprintf('Boundary condition not enforced!\n');
%                     RHO(end,i,3) = ones(size(RHO(end,i,3)));
%                 end
% 
%                 % Total energy equation with unsteady term
%                 PHI(:,i,3) = ((1/dt).*PHI(:,i,2) + eps_2.*phi_laplace(:,i) - (RHO(:,i,2).^(gam-1) - 1)./((gam-1)*M0^2) - 0.5.*(q2_ij(:,i)-1))./((1/(dt)));
%                 PHI(end,i,3) = XX(end,i);
%                 
%             end
%             
%             % Checks
%             if any(RHO(:,i,3) < 0)
%                fprintf('Check density!\nIssue at x = %0.5f\n\n', XX(1,i));
%                fprintf('Setting Negative Density to 0!\n');
%                RHO(RHO(:,i,3) < 0,i,3) = 0;
%             end
%             
%             if (any(abs(phiX(:,i)) > 2)) || (any(abs(phiY(:,i)) > 2))
%                fprintf('Velocity too high!\nIssue at x = %0.5f\n', XX(1,i));
%                fprintf('Max Velocity: %0.5f\n\n', max(max(abs(phiX(:,i))), max(abs(phiY(:,i)))));
%             end
%             
%             if ~isreal(PHI(:,i,3))
%                fprintf('Check phi!\n'); 
%             end
%         end
