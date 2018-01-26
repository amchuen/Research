clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dT = 0.025*pi;
dr = 0.05;

% Field Axis Values
R_max = 15;
r_cyl = 0.5;
T_max = pi;
T_min = 0.0; %(-19*dx);
T_vals = T_min:dT:T_max;
R_vals = r_cyl:dr:R_max;
[TT, RR] = meshgrid(T_vals, R_vals);

XX = RR.*cos(TT);
YY = RR.*sin(TT);

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.45;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_1 = 0.025;
eps_2 = eps_1;
tol = 1e-5;
mintol = 5*tol;
dt = 5e-5;
t_rho = 1;
CFL_on = 1;
iter_min = 800;
show_mesh = 1;

%% BOUNDARY CONDITIONS
BC.Vr_II = cos(TT(end,:)).*(1 - (r_cyl^2)./(RR(end,:).^2));
BC.PHI_II = (RR(end,:) + (r_cyl^2)./(RR(end,:))).*cos(TT(end,:));

%% CALCULATE SOLUTION

PHI = zeros([size(XX), 3]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
PHI(:,:,3) = PHI(:,:,2);
RHO = ones([size(XX), 3]);
res = 1;
ind = [0.0, 0.0];

while (res(end) > tol)|| (length(res) < iter_min) % iterate through time
        
    % Reorganize three-level scheme
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI(:,:,3);
    
    RHO(:,:,1) = RHO(:,:,2);
    RHO(:,:,2) = RHO(:,:,3);
    
    if M0 ~= 0
        % Density Averages
        rho_e = [0.5.*(RHO(:,2:end,2) + RHO(:,1:(end-1),2)), 0.5.*(RHO(:,end,2) + RHO(:,end-1,2))];
        rho_w = [0.5.*(RHO(:,2,2) + RHO(:,1,2)), 0.5.*(RHO(:,1:(end-1),2) + RHO(:,2:end,2))];
        rho_n = [0.5.*(RHO(2:end,:,2) + RHO(1:(end-1),:,2)); 0.5.*(ones(size(RHO(end,:,2))) + RHO(end,:,2))];
        rho_s = 0.5.*([2.*RHO(1,:,2); RHO(2:end,:,2) + RHO(1:(end-1),:,2)]);
        
        % Raidus Averagess
        rr_n = [0.5.*(RR(2:end,:) + RR(1:(end-1),:)); 0.5.*((RR(end,:)+dr) + RR(end,:))];
        rr_s = 0.5.*([2.*RR(1,:); RR(2:end,:) + RR(1:(end-1),:)]);
        
        % Density Gradient -> first order, enforce Neumann conditions at
        % east and west boundaries
        rhoT_w = [(RHO(:,1,2) - RHO(:,2,2))./dT, (RHO(:, 2:end, 2) - RHO(:,1:(end-1),2))./dT];
        rhoT_e = [(RHO(:, 2:end, 2) - RHO(:,1:(end-1),2))./dT, (RHO(:,end-1,2) - RHO(:,end,2))./dT];
        rhoR_n = [(diff(RHO(:,:,2),1,1)./dr); (ones(size(RHO(end,:,2))) - RHO(end,:,2))./dr];
        rhoR_s = [zeros(size(RHO(1,:,2))); (diff(RHO(:,:,2),1,1)./dr)]; % enforce irrotationality at the body
        
        % Velocity -> first order, enforce Neumann conditions at east and
        % west boundaries
        phiT_w = [(PHI(:,1,2) - PHI(:,2,2))./dT, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT];
        phiT_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT, (PHI(:,end-1,2) - PHI(:,end,2))./dT];
        phiR_n = [(diff(PHI(:,:,2),1,1)./dr); BC.Vr_II];
        phiR_s = [zeros(size(PHI(1,:, 2))); (diff(PHI(:,:,2),1,1)./dr)]; % enforce irrotationality at the body

        phiT = 0.5.*(phiT_w + phiT_e);
        phiR = 0.5.*(phiR_n + phiR_s);
        q2_ij = (phiT./RR).^2 + phiR.^2;

        % Mass Flux, Steady Term
        mflux_ss = (rho_e.*phiT_e - rho_w.*phiT_w)./(dT.*RR.^2) + (rr_n.*rho_n.*phiR_n - rr_s.*rho_s.*phiR_s)./(dr.*RR);

        % Laplace of PHI - need to fix here onward
        phi_TT = (phiT_e - phiT_w)./dT + 2.*PHI(:,:,2)./(dT^2);
        phi_RR = (rr_n.*phiR_n - rr_s.*phiR_s)./(dr.*RR) + (rr_n + rr_s).*PHI(:,:,2)./(dr^2);
        
         % Laplace of RHO
        rho_TT = (rhoT_e - rhoT_w)./dT + 2.*RHO(:,:,2)./(dT^2);
        rho_RR = (rr_n.*rhoR_n - rr_s.*rhoR_s)./dr + (rr_n + rr_s).*RHO(:,:,2)./(dr^2);

        % Check CFL conditions
        CFL_i = (max(abs(phiT(:)./(dT.*RR(:).^2))) + max(abs(phiR(:)))./dr)*dt;

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

        % three level scheme activated
        % Total Continuity Equation with unsteady term
        RHO(:,:,3) = -(mflux_ss + (eps_1.*(1./(dT^2.*RR.^2) + 0.5.*(rr_n + rr_s)./(RR.*dr^2)) - 0.5/dt).*RHO(:,:,1) - eps_1.*(rho_TT./(RR.^2) + rho_RR./RR))./(eps_1.*(1./(dT^2.*RR.^2) + 0.5.*(rr_n + rr_s)./(RR.*dr^2)) + 0.5/dt);
        if any(abs((RHO(end,:,3)-1)) > 1e-2)
%             fprintf('Boundary condition not enforced for density!\n');
            RHO(end,abs((RHO(end,:,3)-1)) > 1e-5,3) = 1;
        end

        % Total energy equation with unsteady term
        PHI(:,:,3) = -(((RHO(:,:,2).^(gam-1) - 1)./((gam-1)*M0^2) + 0.5.*(q2_ij - 1)) + (eps_2.*(1./(dT^2.*RR.^2) + 0.5.*(rr_n + rr_s)./(RR.*dr^2)) - 0.5/dt).*PHI(:,:,1) - eps_2.*(phi_TT./(RR.^2) + phi_RR./RR))./(eps_2.*(1./(dT^2.*RR.^2) + 0.5.*(rr_n + rr_s)./(RR.*dr^2)) + 0.5./dt);
        if any(abs((PHI(end,:,3)-BC.PHI_II))./BC.PHI_II > 1e-2)
%             fprintf('Boundary condition not enforced for potential!\n');
%             PHI(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5,3) = XX(end,abs((PHI(end,:,3)-XX(end,:))) > 1e-5);
            PHI(end,:,3) = BC.PHI_II;
        end

        % Checks
        if any(RHO(:,:,3) < 0)
            [r_ind, c_ind] = find(RHO(:,:,3) <0);
           fprintf('Check density!\n');
           fprintf('Issue at %0.2f, %0.2f\n', XX(r_ind,c_ind), YY(r_ind,c_ind));
           fprintf('Setting Negative Density to 0!\n');
           RHO(r_ind, c_ind,3) = 0;
        end

%         if (any(abs(phiT(:)./RR(:)) > 3)) || (any(abs(phiR(:)) > 3))
%             [r_X, c_X] = find(abs(phiT./RR) > 2);
%             [r_Y, c_Y] = find(abs(phiR) > 2);
%            fprintf('Velocity too high!\n');
%            fprintf('Issue at %0.2f, %0.2f\n', XX([r_X;r_Y],[c_X;c_Y]), YY([r_X;r_Y],[c_X; c_Y]));
%            fprintf('Max Velocity: %0.5f\n\n', max(max(abs(phiT)), max(abs(phiR))));
%         end

        if ~isreal(PHI(:,:,3))
           fprintf('Check phi!\n'); 
        end
    end
    
    % Run Residual Check
%     difference = abs((PHI(:,:,3) - PHI(:,:,2))./PHI(:,:,2));
    difference = abs((PHI(:,:,3) - PHI(:,:,2)));
    res(end+1)= max(difference(:));
    [ind(end+1,1), ind(end+1,2)] = find(difference == max(difference(:)), 1);
    
    if (res(end) > 1e1) && (res(end-1) > 1e1)
        figure(1);semilogy(1:length(res), res);
        title(['Residuals, M=' num2str(M0)]);
        xlabel('# of iterations');
        ylabel('Residual (Error)');
        error('Residual growing uncontrollably! Current residual at %0.5f', res(end)); 
    end

    % Initialize next time step
    if (length(res) > 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
    end
    
    % Break if residual starts to grow
    if ((res(end) > res(end-1)) && (res(end) < (mintol))) && (length(res) > iter_min)
        fprintf('Iteration loop breaking at %i iterations.\n', length(res)) 
        break;
    end
    
end

fprintf('\nCalculation Complete!\n');

%% POST-PROCESS

% Calculate Velocities
% phiX = [ones(size(PHI(:,1,3))), (PHI(:,3:end,3) - PHI(:,1:(end-2),3))./(2*dx), ones(size(PHI(:,1,3)))];
% phiY = [dyBdx; (PHI(3:end,:,3) - PHI(1:(end-2),:,3))./(2*dy); zeros(size(dyBdx))];
% q2_ij = phiX.^2 + phiY.^2;

a2_ij = (1./M0.^2)+0.5*(gam-1).*(1 - q2_ij);
M_ij = sqrt(q2_ij./a2_ij);

folderName = ['M_' num2str(M0)];

if ~exist([pwd '\euler_polar\' folderName], 'dir')
    mkdir([pwd '\euler_polar\' folderName]);
end

close all;

% Show Grids
if show_mesh
    figure();
    plot(XX,YY, 'b-', XX', YY', 'b-');
end

figure(1);semilogy(1:length(res), res);
title(['Residual Plot, M=' num2str(M0)]);
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd '\euler_polar\' folderName '\residual_plot.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,round(RHO(:,:,2),3), 50)
title(['Density (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_polar\' folderName '\density.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-(round(q2_ij,4)), 50); %./((RR.*cos(TT)).^2)
title(['Pressure Coefficient Contours, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_polar\' folderName '\cp_contour.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\cp_contour']);

figure();
% plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
plot([fliplr(XX(:,end)'), fliplr(XX(1,:)), XX(:,1)'], 1-[fliplr(q2_ij(:,end)'), fliplr(q2_ij(1,:)), q2_ij(:,1)']);
xlabel('X');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title(['C_p on surface of Cylinder,, M=' num2str(M0)]);
set(gca, 'Ydir', 'reverse');
saveas(gcf, [pwd '\euler_polar\' folderName '\cp_surf.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\cp_surf']);

figure(); % field potential
contourf(XX, YY, round(PHI(:,:,2),5), 50);
title(['Field Potential, \Phi, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_pot.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, round(phiT./RR,5), 50); %./((RR.*cos(TT)).^2)
title(['\Phi_\theta velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_theta.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(XX, YY, phiR, 50); %./((RR.*cos(TT)).^2)
title(['\Phi_R velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_radius.png']);
saveas(gcf, [pwd '\euler_polar\' folderName '\phi_radius']);

% figure();
% contourf(XX, YY, M_ij, 50);
% title(['Mach Number, M=' num2str(M0)]);
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\euler_polar\' folderName '\mach.png']);
% saveas(gcf, [pwd '\euler_polar\' folderName '\mach']);

%% Save Results
save([pwd '\euler_polar\' folderName '\results.mat']);

