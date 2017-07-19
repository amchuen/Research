clc;
clear;
close all;

%% SIM CONTROL PARAMS - FIELD VARIABLE INITIALIZATION 

% Control Params
res = 1;
test_diff = [];
ind = 0;
tol = 0.5e-6;
iter_min = 50000;
enforce_phi = 1;
solve_half = 1;
include_corr = 1;

% Fluid Params
gam = 1.4; % heat 
M0 = 1.1;
v_coeff = 1.25;

% Step Sizes
dT = 0.05*pi;
dr = 0.08;
dt = 1e-5;

alpha = 30;

% Cylinder Dimensions
r_cyl = 0.5;
r_max = 25;

% Field Axis Values
r_vals = (r_cyl+0.5*dr):dr:(r_max);
if M0 > 1 && solve_half
    T_vals = (0.5*pi):dT:pi;
else
    T_vals = 0:dT:(pi);
end

% Grid Setup
[TT, RR] = meshgrid(T_vals, r_vals);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% Field Setup
PHI = ones([size(RR), 2]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
PHI(:,:,3) = PHI(:,:,2);
[n_r, n_T] = size(RR);


BC.Vr_II = cos(TT(end,:)).*(1 - (r_cyl^2)./((RR(end,:) + 0.5*dr).^2));
% BC.PHI_II = (RR(end,:) + (r_cyl^2)./(RR(end,:))).*cos(TT(end,:));
BC.PHI_II = XX(end,:);
visc_on = 0;

%% START SOLUTION
tic
% Raidus Averagess
rr_n = [0.5.*(RR(2:end,:) + RR(1:(end-1),:)); 0.5.*((RR(end,:)+dr) + RR(end,:))];
rr_s = 0.5.*([2.*RR(1,:); RR(2:end,:) + RR(1:(end-1),:)]);

while (res(end) > tol)|| (length(res) < iter_min) % iterate through time
    
    % Initialize next time step
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI(:,:,3);
    
    % Velocity -> first order, enforce Neumann conditions at east and
    % west boundaries
    if M0 > 1 && solve_half
        % outflow boundary condition
        phiT_w = [(PHI(:,2,2) - PHI(:,1,2))./dT, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT];
    else
        % symmetry boundary condition
        phiT_w = [(PHI(:,1,2) - PHI(:,2,2))./dT, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT];
    end
    phiT_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT, (PHI(:,end-1,2) - PHI(:,end,2))./dT]; % symmetry BC
    phiR_n = [(diff(PHI(:,:,2),1,1)./dr); BC.Vr_II];
    phiR_s = [zeros(size(PHI(1,:, 2))); (diff(PHI(:,:,2),1,1)./dr)]; % enforce irrotationality at the body

    phiT = 0.5.*(phiT_w + phiT_e);
    phiR = 0.5.*(phiR_n + phiR_s);
    
    q2_ij = (phiT./RR).^2 + phiR.^2;
	a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
	eps_ij = v_coeff.*max(zeros(size(a2_ij)), 1 - 0.9.*(a2_ij./q2_ij));

    if any(eps_ij(:) ~= 0) && (visc_on ==0) % viscosity not being used and then turned on
        visc_on = 1;
        fprintf('Viscosity model activated! Iteration: %i\n', length(res)+1);        
    elseif (visc_on == 1) && (all(eps_ij(:) == 0))
        visc_on = 0;
        fprintf('Viscosity model turned off! Iteration: %i\n', length(res)+1);
    end
    
%     Check CFL condition
    CFL_i = (max(abs(phiT(:)./(dT.*RR(:)))) + max(abs(phiR(:)))./dr)*dt;
    
    if CFL_i > 1.0
       fprintf('CFL condition not met! Decreasing time steps!\n');
       dt = dt / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
    end
	
	% Initialize Density In-Between Grid
    phi_t = (PHI(:,:,2) - PHI(:,:,1))./dt;
    phi_t_avg = [0.5.*(phi_t(1, 1:(end-1))+ phi_t(1, 2:end));...
				0.25.*(phi_t(2:end, 1:(end-1)) + phi_t(1:(end-1), 1:(end-1)) + phi_t(1:(end-1), 2:end) + phi_t(2:end, 2:end))];
    phiT_avg = [phiT_e(1,1:(end-1)); 0.5.*(phiT_e(2:end,1:(end-1)) + phiT_e(1:(end-1),1:(end-1)))]; % average along radius
    phiR_avg = 0.5*(phiR_s(:,2:end) + phiR_s(:,1:(end-1))); % average along theta direction
	q2_avg = (phiT_avg./rr_s(:,1:(end-1))).^2 + phiR_avg.^2; % total velocity
    a2_avg = (1./M0.^2)-0.5*(gam-1).*(q2_avg - 1); % local speed of sound
	
	rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(q2_avg + 2.*phi_t_avg - 1)).^(1./(gam-1));
    if any(imag(rho_ij(:)) ~= 0)
%         dbstop;
        test = imag(rho_ij);
        [y_ind, x_ind] = find(test ~= 0);
        fprintf('Density is complex at iteration %i!\n', length(res)+1);
        fprintf('Non-real result for density in %i nodes!\n', length(y_ind));
        fprintf('Setting densities to zero!\n\n');
        rho_ij(y_ind, x_ind) = 0;
    end
	
	% Calculate Viscous Dissipation
    del_rho_dT = [0.5.*(rho_ij(:,2)-rho_ij(:,1))-0.5.*(sign(phiT_avg(:,1))).*(rho_ij(:,2)-rho_ij(:,1)),... % Neumann BC b/c of symmetry or outflow at x-axis
            0.5.*(rho_ij(:,3:end)-rho_ij(:,1:(end-2)))-0.5.*(sign(phiT_avg(:,2:(end-1)))).*(rho_ij(:,3:end)-2.*rho_ij(:,2:(end-1))+rho_ij(:,1:(end-2))),... % Neumann BC
            0.5.*(rho_ij(:,end)-rho_ij(:,end-1))-0.5.*(sign(phiT_avg(:,end))).*(-rho_ij(:,end)+rho_ij(:,end-1))];
	
    del_rho_dr = [0.5*(rho_ij(2,:) - rho_ij(1,:))-0.5.*(sign(phiR_avg(1,:))).*(rho_ij(2,:)-rho_ij(1,:));... % Neumann BC due to body
                0.5*(rho_ij(3:end,:) - rho_ij(1:(end-2),:))-0.5.*(sign(phiR_avg(2:(end-1),:))).*(rho_ij(3:end,:)-2.*rho_ij(2:(end-1),:)-rho_ij(1:(end-2),:));...
                0.5*(ones(size(rho_ij(end,:))) - rho_ij(end-1,:))-0.5.*(sign(phiR_avg(end,:))).*(ones(size(rho_ij(end,:)))-2.*rho_ij(end,:)-rho_ij(end-1,:))]; % assumes density to be one at far-field
	rho_ds = ((phiT_avg./rr_n(:,1:(end-1)))./sqrt(q2_avg)).*del_rho_dT + (phiR_avg./sqrt(q2_avg)).*del_rho_dr;
	
	% Calculate Average Densities and Dissipations -> place them on grid-lines
	% theta avg'ed densities should be same size as in-between grid matrix
	RHO_Tavg = [0.5.*(rho_ij(1:(end-1),:) + rho_ij(2:end,:));... % loop through radius, averaged density in radial direction
				0.5.*(rho_ij(end,:) + ones(size(rho_ij(end,:))))];
	RHO_dsT = [0.5.*(rho_ds(1:(end-1),:) + rho_ds(2:end,:));... % loop through radius, averaged density in radial direction
				0.5.*(rho_ds(end,:) + rho_ds(end,:))]; % fixed far-field free-stream condition
%     RHO_dsT = [0.5.*(rho_ds(1:(end-1),:) + rho_ds(2:end,:));... % loop through radius, averaged density in radial direction
% 				0.5.*(rho_ds(end,:) + ones(size(rho_ds(end,:))))]; % Dirichlet condition, fixed far-field free-stream condition
	 
	% radius-avg'ed densities should be same size as original grid
	RHO_Ravg = [0.5.*(rho_ij(:,1) + rho_ij(:,1)),... % Neumann condition at TE line, may need to change if 
				0.5.*(rho_ij(:,2:(end)) + rho_ij(:,1:(end-1))),...
				0.5.*(rho_ij(:,end) + rho_ij(:,end))]; % Neumann condition at LE line
	RHO_dsR = [0.5.*(rho_ds(:,1) + rho_ds(:,1)),... % Neumann condition at TE line, may need to change if 
				0.5.*(rho_ds(:,2:(end)) + rho_ds(:,1:(end-1))),...
				0.5.*(rho_ds(:,end) + rho_ds(:,end))]; % Neumann condition at LE line
            
    % Apply Elliptic Stencil with Three-Level Scheme
    RHO_W = [RHO_Tavg(:,1),RHO_Tavg] - eps_ij.*[RHO_dsT(:,1),RHO_dsT];
    RHO_E = [RHO_Tavg,RHO_Tavg(:,end)] - eps_ij.*[RHO_dsT,RHO_dsT(:,end)];
    RHO_N = [RHO_Ravg(2:end, :) - eps_ij(1:(end-1),:).*RHO_dsR(2:end, :);...
            ones(size(RHO_Ravg(end,:)))];
    RHO_S = RHO_Ravg - eps_ij.*RHO_dsR;
    
    PHI_RR = (rr_n.*RHO_N.*phiR_n - rr_s.*RHO_S.*phiR_s)./dr;
    PHI_TT = (RHO_E.*phiT_e - RHO_W.*phiT_w)./dT;
%     RHO = (1 - 0.5*(gam-1).*M0.^2.*(q2_ij - 1)).^(1./(gam-1));
    if include_corr
        RHO = 0.25.*([RHO_Ravg(2:end, :); ones(size(RHO_Ravg(end,:)))] + ...
                     RHO_Ravg + [RHO_Tavg(:,1),RHO_Tavg] + ...
                     [RHO_Tavg,RHO_Tavg(:,end)]);
    %     RHO = 0.25.*(RHO_W + RHO_E + RHO_N + RHO_S);
        CORR = ((RHO.^(gam - 1) - 1)./((gam-1)*M0^2) + 0.5.*(q2_ij - 1));
    else
        CORR = zeros(size(RHO_W));
    end

    PHI(:,:,3) = ((2./(dt^2) - alpha/dt).*PHI(:,:,2) - (1/dt^2 - alpha/dt).*PHI(:,:,1) + PHI_RR./RR + PHI_TT./(RR.^2) - alpha.*CORR)./(1/dt^2);
    if enforce_phi
%         PHI(end,:,3) = XX(end,:);
        PHI(end,:,3) = BC.PHI_II;
    end

    % Run Residual Check
    difference = abs(PHI(:,:,3) - PHI(:,:,2));
    [res(end+1), ind(end+1)] = max(difference(:));
    
    if (length(res) > 500) && (mod(length(res), 2000) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
        toc;
        figure(1);semilogy(1:length(res), res);
        fprintf('\n');
    end
    
    if ~isreal(PHI(:,:,end)) || any(any(isnan(PHI(:,:,end))))
        fprintf('PHI exhibits non-solutions (either non-real or NaN) in nodes!\n');
    end

end

toc
%% Post Process
if M0 > 1 && solve_half
    % outflow boundary condition
    phiT_w = [(PHI(:,2,2) - PHI(:,1,2))./dT, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT];
else
    % symmetry boundary condition
    phiT_w = [(PHI(:,1,2) - PHI(:,2,2))./dT, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT];
end
phiT_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./dT, (PHI(:,end-1,2) - PHI(:,end,2))./dT]; % symmetry BC
phiR_n = [(diff(PHI(:,:,2),1,1)./dr); BC.Vr_II];
phiR_s = [zeros(size(PHI(1,:, 2))); (diff(PHI(:,:,2),1,1)./dr)]; % enforce irrotationality at the body

phiT = 0.5.*(phiT_w + phiT_e);
phiR = 0.5.*(phiR_n + phiR_s);

q2_ij = (phiT./RR).^2 + phiR.^2;
a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
M2_ij = q2_ij./a2_ij; % local Mach number
M_ij = sqrt(M2_ij);

folderName = ['M_' num2str(M0)];
type = '\cyl_mod\' ;

if ~exist([pwd type folderName], 'dir')
    mkdir([pwd type folderName]);
end

close all;
figure(1);semilogy(1:length(res), res);
title('Residual Plot');
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd type folderName '\residual_plot.png']);
saveas(gcf, [pwd type folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,RHO_Ravg, 50)
title('Density (Normalized)');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd type folderName '\density.png']);
saveas(gcf, [pwd type folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-(q2_ij), 50); %./((RR.*cos(TT)).^2)
% set(hC,'LineStyle','none');

% hold on;
% [hC hC] = contourf(XX, -1.*YY, 1-q2_ij, 250);
% set(hC,'LineStyle','none');

title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd type folderName '\cp_contour.png']);
saveas(gcf, [pwd type folderName '\cp_contour']);

% figure();
% % plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
% plot(fliplr(TT(1,:)).*180 ./ pi, 1-q2_ij(1,:));
% xlabel('\theta');
% %     ylabel('\phi_{\theta}');
% ylabel('C_p');
% title('C_p on surface of Cylinder');
% % set(gca, 'Xdir', 'reverse');
% saveas(gcf, [pwd type folderName '\cp_surf.png']);
% saveas(gcf, [pwd type folderName '\cp_surf']);

figure();
if M0>1 && solve_half
    plot([fliplr(XX(:,end)'), fliplr(XX(1,:))], 1-[fliplr(q2_ij(:,end)'), fliplr(q2_ij(1,:))]); 
else
    plot([fliplr(XX(:,end)'), fliplr(XX(1,:)), XX(:,1)'], 1-[fliplr(q2_ij(:,end)'), fliplr(q2_ij(1,:)), q2_ij(:,1)']);
end

% plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);

xlabel('X');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title(['C_p on surface of Cylinder,, M=' num2str(M0)]);
set(gca, 'Ydir', 'reverse');
% axis equal
saveas(gcf, [pwd type folderName '\cp_surf.png']);
saveas(gcf, [pwd type folderName '\cp_surf']);

figure(); % field potential
contourf(XX, YY, PHI(:,:,end), 50);
title('Field Potential, \Phi');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd type folderName '\phi_pot.png']);
saveas(gcf, [pwd type folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, phiT./RR, 50); %./((RR.*cos(TT)).^2)
title('\Phi_\theta velocity');
colorbar('eastoutside');
saveas(gcf, [pwd type folderName '\phi_theta.png']);
saveas(gcf, [pwd type folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(XX, YY, phiR, 50); %./((RR.*cos(TT)).^2)
title('\Phi_R velocity');
colorbar('eastoutside');
saveas(gcf, [pwd type folderName '\phi_radius.png']);
saveas(gcf, [pwd type folderName '\phi_radius']);

figure();
contourf(XX, YY, sqrt(M2_ij), 50);
title('Mach Number');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd type folderName '\mach.png']);
saveas(gcf, [pwd type folderName '\mach']);

%% Save Results
save([pwd type folderName '\results.mat']);

%% Old Code
% 	Apply Elliptic Stencil
%     for i = 1:length(T_vals) % we will not calculate the far-field b/c that is already defined
%         % can be solved marching w/ flow or against flow
%         % currently set up to march against flow, if marching with flow
%         % west becomes east, and east becomes west
% 
%         % Calculate Radius Averages
%         RR_avg = 0.5.*(RR(2:end,i) + RR(1:(end-1),i));
% 
%         % Calculate Density Averages
%         if i == 1
%             RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
%             RHO_E = RHO_W;
%         elseif i == length(T_vals)
%             RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
%             RHO_W = RHO_E;
%         else
%             RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
%             RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
%         end
% 
%         RHO_N = RHO_Ravg(2:end, i) - v_coeff.*eps_ij(1:(end-1),i).*RHO_dsR(2:end, i);
%         RHO_S = RHO_Ravg(1:(end-1), i) - v_coeff.*eps_ij(1:(end-1),i).*RHO_dsR(1:(end-1), i);
% 
%         PHI_RR = [(RR_avg(1)*RHO_N(1)*(PHI(2,i,2) - PHI(1,i,2))/dr)/(RR(1,i)*dr);...
%                 (RR_avg(2:end).*RHO_N(2:end).*(PHI(3:end,i,2) - PHI(2:(end-1),i,2))./dr - RR_avg(1:(end-1)).*RHO_S(2:end).*(PHI(2:(end-1),i,2) - PHI(1:(end-2),i,2))./dr)./(RR(2:(end-1),i).*dr)];
% %         PHI_RR = [(rr_n(1,i).*RHO_N(1)*phiR_n(1,i) - rr_s(1,i).*RHO_S(1).*phiR_s(1,i));...
% %                     (rr_n(2:(end-1),i).*RHO_N(2:end)*phiR_n(2:(end-1),i) - rr_s(2:(end-1),i).*RHO_S(2:end).*phiR_s(2:(end-1),i))]./(RR(2:(end-1),i)*dr);
% 
%         if i == 1
%             PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i+1,2))./dT);
%         elseif i == length(T_vals)
%             PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i-1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
%         else
%             PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
%         end
% 
%         PHI(1:(end-1), i, 3) = (PHI_RR + PHI_TT - alpha.*CORR(1:(end-1),i) + 2./(dt^2).*PHI(1:(end-1),i,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,2-1))./(1/dt^2 + 0.5*alpha/dt);
% 
%     end