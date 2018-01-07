clc;
clear;
close all;

%% SIM CONTROL PARAMS - FIELD VARIABLE INITIALIZATION 

% Fluid Params
gam = 1.4; % heat 
M0 = 0.3;
v_coeff = 1.5;

% Step Sizes
dr = 0.1;
dT = 0.005*pi;
dt = 0.01;

alpha = 30;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 8.0;

% Field Axis Values
r_vals = (r_cyl+0.5*dr):dr:(r_max);
if M0 >=1.0
    T_vals = (0.5*pi):dT:(pi);  
else
%     T_vals = 0:dT:(pi);  
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
res = 1;
ind = 0;
tol = 1e-5;

BC.Vr_II = cos(TT(end,:)).*(1 - (r_cyl^2)./(RR(end,:).^2));
BC.dirichlet.PHI_II = (RR(end,:) + (r_cyl^2)./(RR(end,:))).*cos(TT(end,:));
visc_on = 0;

%% START SOLUTION
tic
while res(end) > tol % iterate through time
    
    % Initialize next time step
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI(:,:,3);
    
    if (length(res) > 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
        toc;
        figure(1);semilogy(1:length(res), res);
        fprintf('\n');
    end

%     PHI(:,:,iter+1) = PHI(:,:,iter); % setup n+1th time step
    
    if ~isreal(PHI(:,:,end)) || any(any(isnan(PHI(:,:,end))))
        fprintf('PHI exhibits non-solutions (either non-real or NaN) in nodes!\n');
    end
    
	% Get Artificial Visc at Nodes
    if M0 >= 1.0
%         phiT = [(PHI(:,2,end) - PHI(:,1,end))./(dT .* RR(:,1)),...
%             (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*dT .* RR(:, 2:(end-1))),... % Central Difference
%             zeros(size(RR(:,end)))]; % PHI_X    % Neumann condition at LE/TE of cyl b/c of symmetry condition

        phiT = [zeros(size(RR(:,1))),... % Neumann condition at LE/TE of cyl b/c of symmetry condition
            (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*dT .* RR(:, 2:(end-1))),... % Central Difference
            (PHI(:,end,end) - PHI(:,end-1,end))./(dT .* RR(:,end))]; % PHI_X
    else
        phiT = [zeros(size(RR(:,1))),... % Neumann condition at LE/TE of cyl b/c of symmetry condition
            (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*dT .* RR(:, 2:(end-1))),... % Central Difference
            zeros(size(RR(:,end)))]; % PHI_X
    end
	phiR  = [zeros(size(PHI(1,:, 2)));... % Neumann condition at cylinder surface
			(PHI(3:end,:,2) - PHI(1:(end-2),:,2))./(2*dr);...
			BC.Vr_II]; % PHI_Y, central difference
	
	q2_ij = phiT.^2 + phiR.^2;
	a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
	eps_ij = max(zeros(size(a2_ij)), 1 - 0.99.*(a2_ij./q2_ij));

    if any(eps_ij(:) ~= 0) && (visc_on ==0) % viscosity not being used and then turned on
        visc_on = 1;
        fprintf('Viscosity model activated! Iteration: %i\n', length(res)+1);        
    elseif (visc_on == 1) && (all(eps_ij(:) == 0))
        visc_on = 0;
        fprintf('Viscosity model turned off! Iteration: %i\n', length(res)+1);
    end
    
    % Check CFL condition
    CFL_i = (max(abs(phiT(:)./(dT.*RR(:)))) + max(abs(phiR(:)))./dr)*dt;
    
%     if CFL_i > 1.0
%        fprintf('CFL condition not met! Decreasing time steps!\n');
%        dt = dt / CFL_i;
%        fprintf('New time step:%0.5f\n', dt);
%     end
	
	% Initialize Uncorrected Density
    phi_t = (PHI(:,:,2) - PHI(:,:,1))./dt;
    phi_t_avg = [0.5.*(phi_t(1, 1:(end-1))+ phi_t(1, 2:end));...
				0.25.*(phi_t(2:end, 1:(end-1)) + phi_t(1:(end-1), 1:(end-1)) + phi_t(1:(end-1), 2:end) + phi_t(2:end, 2:end))];
	u_avg = [(PHI(1,2:end,2) - PHI(1,1:(end-1),2))./(RR(1,1:(end-1)).*dT);...
			0.5.*((PHI(2:end,2:end,2) - PHI(2:end,1:(end-1),2))./(RR(2:end,1:(end-1)).*dT) + (PHI(1:(end-1),2:(end),2) - PHI(1:(end-1),1:(end-1),2))./(RR(1:(end-1),1:(end-1)).*dT))];
	v_avg = [zeros(1,n_T-1);...
			0.5.*((PHI(2:end,1:(end-1),2) - PHI(1:(end-1),1:(end-1),2))./dr + (PHI(2:end,2:end,2) - PHI(1:(end-1),2:end,2))./dr)];
	
	q2_avg = u_avg.^2 + v_avg.^2; % total velocity
    a2_avg = (1./M0.^2)-0.5*(gam-1).*(q2_avg - 1);
	
	rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(q2_avg + 2.*phi_t_avg - 1)).^(1./(gam-1));
    if ~isreal(rho_ij)
        fprintf('Density is complex at iteration %i!\n', length(res)+1);
    end
	
	% Calculate Viscous Dissipation
	del_rho_dT = [0.5.*(rho_ij(:,2)-rho_ij(:,1))-0.5.*(sign(u_avg(:,1))).*(rho_ij(:,2)-rho_ij(:,1)),... % Neumann BC b/c of symmetry at x-axis
			0.5.*(rho_ij(:,3:end)-rho_ij(:,1:(end-2)))-0.5.*(sign(u_avg(:,2:(end-1)))).*(rho_ij(:,3:end)-2.*rho_ij(:,2:(end-1))+rho_ij(:,1:(end-2))),... % Neumann BC
			0.5.*(rho_ij(:,end)-rho_ij(:,end-1))-0.5.*(sign(u_avg(:,end))).*(-rho_ij(:,end)+rho_ij(:,end-1))];
	
	del_rho_dr = [0.5*(rho_ij(2,:) - rho_ij(1,:))-0.5.*(sign(v_avg(1,:))).*(rho_ij(2,:)-rho_ij(1,:));... % Neumann BC due to body
			0.5*(rho_ij(3:end,:) - rho_ij(1:(end-2),:))-0.5.*(sign(v_avg(2:(end-1),:))).*(rho_ij(3:end,:)-2.*rho_ij(2:(end-1),:)-rho_ij(1:(end-2),:));...
			0.5*(ones(size(rho_ij(end,:))) - rho_ij(end-1,:))-0.5.*(sign(v_avg(end,:))).*(ones(size(rho_ij(end,:)))-2.*rho_ij(end,:)-rho_ij(end-1,:))]; % assumes density to be one at far-field
	
	rho_ds = (u_avg./sqrt(q2_avg)).*del_rho_dT + (v_avg./sqrt(q2_avg)).*del_rho_dr;
	
	% Calculate Average Densities and Dissipations -> place them on grid-lines
	% theta avg'ed densities should be same size as in-between grid matrix
	RHO_Tavg = [0.5.*(rho_ij(1:(end-1),:) + rho_ij(2:end,:));... % loop through radius, averaged density in radial direction
				0.5.*(rho_ij(end,:) + ones(size(rho_ij(end,:))))];
	RHO_dsT = [0.5.*(rho_ds(1:(end-1),:) + rho_ds(2:end,:));... % loop through radius, averaged density in radial direction
				0.5.*(rho_ds(end,:) + ones(size(rho_ds(end,:))))]; % Dirichlet condition, fixed far-field free-stream condition
	 
	% radius-avg'ed densities should be same size as original grid
	RHO_Ravg = [0.5.*(rho_ij(:,1) + rho_ij(:,1)),... % Neumann condition at TE line, may need to change if 
				0.5.*(rho_ij(:,2:(end)) + rho_ij(:,1:(end-1))),...
				0.5.*(rho_ij(:,end) + rho_ij(:,end))]; % Neumann condition at LE line
	RHO_dsR = [0.5.*(rho_ds(:,1) + rho_ds(:,1)),... % Neumann condition at TE line, may need to change if 
				0.5.*(rho_ds(:,2:(end)) + rho_ds(:,1:(end-1))),...
				0.5.*(rho_ds(:,end) + rho_ds(:,end))]; % Neumann condition at LE line
				
	% Apply Elliptic Stencil
    if M0 >=1.0
        % for the supersonic case, we will be expecting a bow-shock...
        % thus we will need to extrapolate the set of points at the end of
        % grid, which is the location above the cylinder
        
        % because the flow is supersonic, we will not need to worry about
        % boundary conditions and transfering of params
        
        for i = 1:length(T_vals) % we will not calculate the far-field b/c that is already defined
%             % Calculate Radius Averages
%             RR_avg = 0.5.*(RR(2:end,i) + RR(1:(end-1),i));
% 
%             % Calculate Density Averages
%             if i == 1 % extrapolate from forward two points
%                 RHO_W = RHO_Tavg(1:(end-1),i+1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i+1);
%                 RHO_E = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
%             elseif i == length(T_vals) % symmetry at leading edge
%                 RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
%                 RHO_W = RHO_E;
%             else
%                 RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
%                 RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
%             end
% 
%             RHO_N = RHO_Ravg(2:end, i) - v_coeff.*eps_ij(1:(end-1),i).*RHO_dsR(2:end, i);
%             RHO_S = RHO_Ravg(2:(end-1), i) - v_coeff.*eps_ij(2:(end-1),i).*RHO_dsR(2:(end-1), i);
% 
%             PHI_RR = [(RR_avg(1)*RHO_N(1)*(PHI(2,i,2) - PHI(1,i,2))/dr)/(RR(1,i)*dr);...
%                     (RR_avg(2:end).*RHO_N(2:end).*(PHI(3:end,i,2) - PHI(2:(end-1),i,2))./dr - RR_avg(1:(end-1)).*RHO_S.*(PHI(2:(end-1),i,2) - PHI(1:(end-2),i,2))./dr)./(RR(2:(end-1),i).*dr)];
% 
%             if i == 1 % symmetry at leading edge
% %                 PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT - RHO_W.*(PHI(1:(end-1),i-1,2) - PHI(1:(end-1),i-2,2))./dT);
%                 PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+2,2) - PHI(1:(end-1),i+1,2))./dT - RHO_E.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT);
%             elseif i == length(T_vals) % extrapolate from previous two points
%                 PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i-1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
%             else
%                 PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
%             end
% 
%             if i == 1%length(T_vals) % enforce velocity propagation at outlet? PHI_TT = 0?
%                 PHI(1:(end-1), i, 3) = ((2/(dt^2)).*(2.*PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i+2,2)) -(1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,1))./(1/dt^2 + 0.5*alpha/dt);
% %                 PHI(1:(end-1), i, 3) = (PHI_TT + 2./(dt^2).*PHI(1:(end-1),i,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,2-1))./(1/dt^2 + 0.5*alpha/dt);
%             else
%                 PHI(1:(end-1), i, 3) = (PHI_RR + PHI_TT + 2./(dt^2).*PHI(1:(end-1),i,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,2-1))./(1/dt^2 + 0.5*alpha/dt);
%             end

            % Calculate Radius Averages
            RR_avg = 0.5.*(RR(2:end,i) + RR(1:(end-1),i));

            % Calculate Density Averages
            if i == 1
                RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
                RHO_E = RHO_W;
            elseif i == length(T_vals)
                RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
                RHO_W = RHO_E;
            else
                RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
                RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
            end

            RHO_N = RHO_Ravg(2:end, i) - v_coeff.*eps_ij(1:(end-1),i).*RHO_dsR(2:end, i);
            RHO_S = RHO_Ravg(2:(end-1), i) - v_coeff.*eps_ij(2:(end-1),i).*RHO_dsR(2:(end-1), i);

            PHI_RR = [(RR_avg(1)*RHO_N(1)*(PHI(2,i,2) - PHI(1,i,2))/dr)/(RR(1,i)*dr);...
                    (RR_avg(2:end).*RHO_N(2:end).*(PHI(3:end,i,2) - PHI(2:(end-1),i,2))./dr - RR_avg(1:(end-1)).*RHO_S.*(PHI(2:(end-1),i,2) - PHI(1:(end-2),i,2))./dr)./(RR(2:(end-1),i).*dr)];

            if i == 1
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i+1,2))./dT);
            elseif i == length(T_vals)
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i-1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
            else
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
            end

            PHI(1:(end-1), i, 3) = (PHI_RR + PHI_TT + 2./(dt^2).*PHI(1:(end-1),i,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,2-1))./(1/dt^2 + 0.5*alpha/dt);

        end
        PHI(end,:, 3) = BC.dirichlet.PHI_II; % apply BC at outer boundary        
        % now extrapolate 1st line from forward two lines using elliptic
        % scheme... do we need to correct with time?
        
%         RHO_W = RHO_Tavg(1:(end-1),2) - e_coeff.*eps_ij(1:(end-1), 1).*RHO_dsT(1:(end-1), 2);
%         RHO_E = RHO_Tavg(1:(end-1),1) - e_coeff.*eps_ij(1:(end-1), 1).*RHO_dsT(1:(end-1), 1);
%         
%         RHO_N = RHO_Ravg(2:end, 2) - e_coeff.*eps_ij(1:(end-1),1).*RHO_dsR(2:end, 2);
%         RHO_S = RHO_Ravg(2:(end-1), 2) - e_coeff.*eps_ij(2:(end-1),1).*RHO_dsR(2:(end-1), 2);
%         
%         PHI_RR = [(RR_avg(1)*RHO_N(1)*(PHI(2,2,2) - PHI(1,2,2))/dr)/(RR(1,2)*dr);...
%                     (RR_avg(2:end).*RHO_N(2:end).*(PHI(3:end,2,2) - PHI(2:(end-1),2,2))./dr - RR_avg(1:(end-1)).*RHO_S.*(PHI(2:(end-1),2,2) - PHI(1:(end-2),2,2))./dr)./(RR(2:(end-1),2).*dr)];
%                 
%         PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
%         
%         PHI(1:(end-1), 1, 3) = (PHI_RR + PHI_TT + 2./(dt^2).*PHI(1:(end-1),1,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),1,2-1))./(1/dt^2 + 0.5*alpha/dt);
        
        
        
    else 
        for i = 1:length(T_vals) % we will not calculate the far-field b/c that is already defined
            % can be solved marching w/ flow or against flow
            % currently set up to march against flow, if marching with flow
            % west becomes east, and east becomes west
            
            % Calculate Radius Averages
            RR_avg = 0.5.*(RR(2:end,i) + RR(1:(end-1),i));

            % Calculate Density Averages
            if i == 1
                RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
                RHO_E = RHO_W;
            elseif i == length(T_vals)
                RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
                RHO_W = RHO_E;
            else
                RHO_W = RHO_Tavg(1:(end-1),i) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i);
                RHO_E = RHO_Tavg(1:(end-1),i-1) - v_coeff.*eps_ij(1:(end-1), i).*RHO_dsT(1:(end-1), i-1);
            end

            RHO_N = RHO_Ravg(2:end, i) - v_coeff.*eps_ij(1:(end-1),i).*RHO_dsR(2:end, i);
            RHO_S = RHO_Ravg(2:(end-1), i) - v_coeff.*eps_ij(2:(end-1),i).*RHO_dsR(2:(end-1), i);

            PHI_RR = [(RR_avg(1)*RHO_N(1)*(PHI(2,i,2) - PHI(1,i,2))/dr)/(RR(1,i)*dr);...
                    (RR_avg(2:end).*RHO_N(2:end).*(PHI(3:end,i,2) - PHI(2:(end-1),i,2))./dr - RR_avg(1:(end-1)).*RHO_S.*(PHI(2:(end-1),i,2) - PHI(1:(end-2),i,2))./dr)./(RR(2:(end-1),i).*dr)];

            if i == 1
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i+1,2))./dT);
            elseif i == length(T_vals)
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i-1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
            else
                PHI_TT = (1./RR(1:(end-1),i).^2).*(1/dT).*(RHO_W.*(PHI(1:(end-1),i+1,2) - PHI(1:(end-1),i,2))./dT - RHO_E.*(PHI(1:(end-1),i,2) - PHI(1:(end-1),i-1,2))./dT);
            end

            PHI(1:(end-1), i, 3) = (PHI_RR + PHI_TT + 2./(dt^2).*PHI(1:(end-1),i,2) - (1/dt^2 - 0.5*alpha/dt).*PHI(1:(end-1),i,2-1))./(1/dt^2 + 0.5*alpha/dt);

        end
        PHI(end,:, 3) = BC.dirichlet.PHI_II; % apply BC at outer boundary
    end
    
    % Run Residual Check
    difference = abs(PHI(:,:,3) - PHI(:,:,2));
    [res(end+1), ind(end+1)] = max(difference(:));
    
end

toc
%% Post Process

% Calculate Velocity as function of PHI(n)
if M0 < 0.8
   phiT = [zeros(size(RR(:,1))),... % Neumann condition at LE/TE of cyl b/c of symmetry condition
        (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*dT .* RR(:, 2:(end-1))),... % Central Difference
        zeros(size(RR(:,end)))]; % PHI_X
else
   phiT = [(PHI(:,2,end) - PHI(:,1,end))./(dT .* RR(:,1)),...
        (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*dT .* RR(:, 2:(end-1))),... % Central Difference
        zeros(size(RR(:,end)))]; % PHI_X    % Neumann condition at LE/TE of cyl b/c of symmetry condition
end
phiR  = [zeros(size(PHI(1,:,end)));... % Neumann condition at cylinder surface
		(PHI(3:end,:,end) - PHI(1:(end-2),:,end))./(2*dr);...
		BC.Vr_II]; % PHI_Y, central difference

q2_ij = phiT.^2 + phiR.^2;
a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
M2_ij = q2_ij./a2_ij; % local Mach number
M_ij = sqrt(M2_ij);

folderName = ['M_' num2str(M0)];

if ~exist([pwd '\cylinder\' folderName], 'dir')
    mkdir([pwd '\cylinder\' folderName]);
end

close all;
figure(1);semilogy(1:length(res), res);
title('Residual Plot');
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd '\cylinder\' folderName '\residual_plot.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,RHO_Ravg, 50)
title('Density (Normalized)');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\cylinder\' folderName '\density.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-(q2_ij), 50); %./((RR.*cos(TT)).^2)
title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\cylinder\' folderName '\cp_contour.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\cp_contour']);

figure();
% plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
plot(fliplr(TT(1,:)).*180 ./ pi, 1-q2_ij(1,:));
xlabel('\theta');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title('C_p on surface of Cylinder');
% set(gca, 'Xdir', 'reverse');
saveas(gcf, [pwd '\cylinder\' folderName '\cp_surf.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\cp_surf']);

figure(); % field potential
contourf(XX, YY, PHI(:,:,end), 50);
title('Field Potential, \Phi');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\cylinder\' folderName '\phi_pot.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, phiT, 50); %./((RR.*cos(TT)).^2)
title('\Phi_\theta velocity');
colorbar('eastoutside');
saveas(gcf, [pwd '\cylinder\' folderName '\phi_theta.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(XX, YY, phiR, 50); %./((RR.*cos(TT)).^2)
title('\Phi_R velocity');
colorbar('eastoutside');
saveas(gcf, [pwd '\cylinder\' folderName '\phi_radius.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\phi_radius']);

figure();
contourf(XX, YY, sqrt(M2_ij), 50);
title('Mach Number');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\cylinder\' folderName '\mach.png']);
saveas(gcf, [pwd '\cylinder\' folderName '\mach']);

%% Save Results
save([pwd '\cylinder\' folderName '\results.mat']);