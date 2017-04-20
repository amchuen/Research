clc;
clear;
close all;

%% SIM CONTROL PARAMS - GRID INITIALIZATION

% Step Sizes
dr = 0.1;
dT = 0.01*pi;
dt = 0.1*dr;

alpha = 2.5;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 10;

% Field Axis Values
r_vals = (r_cyl+0.5*dr):dr:(r_max);
T_vals = 0:dT:(pi);

% Time Inputs
% start = 0.0;
% stop = 100 * run.times.dt;
% dt = [0.001];

% Grid Setup
[TT, RR] = meshgrid(T_vals, r_vals);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% Fluid Params
gam = 1.4; % heat 
M0 = 0.25;
visc = 0.0;

%% FIELD VARIABLE INITIALIZATION

PHI = ones([size(RR), 2]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
[n_r, n_T] = size(RR);
res = 1;
ind = 0;
iter = 1.0;
tol = 1e-5;

BC.dirichlet.Vr_II = cos(TT(end,:)).*(1 - (r_cyl^2)./(RR(end,:).^2));
BC.dirichlet.PHI_II = (RR(end,:) + (r_cyl^2)./(RR(end,:))).*cos(TT(end,:));

%% START SOLUTION

while res(end) > tol % iterate through time
    
    % Initialize next time step
    if (iter > 500) && (mod(iter, 500) == 0)
        fprintf('Iteration Ct: %i\n', iter);
        fprintf('Current Residual: %0.5f\n', res(end));
        figure(1);semilogy(1:iter, res);
    end
    iter = iter + 1; % use this to call nth time step
    PHI(:,:,iter+1) = PHI(:,:,iter); % setup n+1th time step
    
    if ~isreal(PHI(:,:,iter)) || any(any(isnan(PHI(:,:,iter))))
        fprintf('PHI exhibits non-solutions (either non-real or NaN) in nodes!\n');
    end

    % Calculate Velocity as function of PHI(n)
    U_n = zeros(size(RR)); % PHI_X
    V_n = zeros(size(RR)); % PHI_Y
    for i = 2:(size(PHI,2)-1) % loop through thetas
        U_n(:,i) = (PHI(:,i+1,iter) - PHI(:,i-1,iter))./(2*dT .* RR(:, i));
    end
    
    V_n(end,:) = BC.dirichlet.Vr_II; % apply far-field potential flow condition
    for ii = 2:(size(PHI,1)-1) % loop through radii
        V_n(ii,:) = (PHI(ii+1,:,iter) - PHI(ii-1,:,iter))./(dr); % central difference
    end
    
    % Get Local Viscosity
    q2_ij = U_n.^2 + V_n.^2;
    a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
    M2_ij = q2_ij./a2_ij; % local Mach number
    
    
    % Initialize Density
    phi_t = (PHI(:,:,iter) - PHI(:,:,iter-1))./dt;
    rho_ij = zeros([n_r, n_T-1]);
    u_avg = rho_ij;
    v_avg = rho_ij;
    phi_t_avg = rho_ij;
    
    for j = 2:(n_r) % assume same number of points in the radius direction
        % looping through radius points
        % here, we want to determine density through all the radii except
        % for at the cylinder first, to account for the offset grid
        
        for ii = 1:(n_T-1) % assume one less set of points in theta direction, due to between-grid calculation
            % looping through angular points
            % Will have 1 less point than grid b/c of offset
            % Neumann enforced for density by setting end grid-points to be
            % same if reflecte across boundary
            
            % Calculate average velocities...
            u_avg(j,ii) = 0.5*((PHI(j,ii+1,iter) - PHI(j,ii,iter))/(RR(j,ii)*dT) + (PHI(j-1,ii+1,iter) - PHI(j-1,ii,iter))/(RR(j-1,ii)*dT));
            v_avg(j,ii) = 0.5*((PHI(j,ii,iter) - PHI(j-1,ii,iter))/dr + (PHI(j,ii+1,iter) - PHI(j-1,ii+1,iter))/dr);
            phi_t_avg(j,ii) = mean([phi_t(j, ii), phi_t(j-1, ii), phi_t(j-1, ii+1), phi_t(j, ii+1)]);            
        end
    end
    
    % Apply Neumann condition on velocities
    u_avg(1,:) = u_avg(2,:);
    v_avg(1,:) = v_avg(2,:);
    phi_t_avg(1,:) = phi_t_avg(2,:);
    
    % Calculate Local (Artificial) Viscosity
    q2_avg = u_avg.^2 + v_avg.^2; % total velocity
    a2_avg = (1./M0.^2)-0.5*(gam-1).*(q2_avg - 1);
    eps_ij = 1.3.*max(zeros(size(rho_ij)), 1 - 0.99.*(a2_avg./q2_avg));
    
    % Calculate uncorrected density
    test =  q2_avg + 2.*phi_t_avg - 1;
    if (1 - 0.5*(gam-1).*M0^2*(max(abs(test(:))))) < 0
        fprintf('Density is complex! Please resolve before continuing!\n');
    end
    
    rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(q2_avg + 2.*phi_t_avg - 1)).^(1./(gam-1));
    
    if (M0 == 0) & ((rho_ij(1,:) ~= 1))
        fprintf('Incompressible Case not satisfied!\n');
    end
    
    % Calculate viscosity corrections
    del_rho_dT = zeros(size(rho_ij));
    for i = 1:(size(rho_ij,2)) % loop through theta
        if i ==1 % Neumann BC
            i1 = i+1;
            i_1 = i;
        elseif i == size(rho_ij,2) % Neumann BC
            i1 = i;
            i_1 = i-1;
        else
            i1 = i+1;
            i_1 = i-1;
        end
        del_rho_dT(:,i) = 0.5.*(rho_ij(:,i1)-rho_ij(:,i_1))-0.5.*(sign(u_avg(:,i))).*(rho_ij(:,i1)-2.*rho_ij(:,i)+rho_ij(:,i_1));
    end
    
    del_rho_dr = zeros(size(rho_ij));
    for ii = 1:(size(rho_ij,1)) % loop through radius
        if ii ==1 % Neumann BC
            del_rho_dr(ii,:) = 0.5*(rho_ij(ii+1,:) - rho_ij(ii,:))-0.5.*(sign(v_avg(ii,:))).*(rho_ij(ii+1,:)-rho_ij(ii,:));
        elseif ii == size(rho_ij,1) % Dirichlet BC
            del_rho_dr(ii,:) = 0.5*(ones(size(rho_ij(ii,:))) - rho_ij(ii,:))-0.5.*(sign(v_avg(ii,:))).*(ones(size(rho_ij(ii,:)))-2.*rho_ij(ii,:)-rho_ij(ii-1,:)); % assumes density to be zero at far-field
        else
            del_rho_dr(ii,:) = 0.5*(rho_ij(ii+1,:) - rho_ij(ii,:))-0.5.*(sign(v_avg(ii,:))).*(rho_ij(ii+1,:)-2.*rho_ij(ii,:)-rho_ij(ii-1,:));
        end
    end
    
    rho_ds = (u_avg./sqrt(q2_avg)).*del_rho_dT + (v_avg./sqrt(q2_avg)).*del_rho_dr;
    RHO = rho_ij - eps_ij.*rho_ds;
    
    % Calculate Average Densities on Grid
    RHO_Tavg = zeros(size(RHO));
    for i = 1:size(RHO,1) % loop through radius
        if i == size(RHO,1)
            RHO_Tavg(i,:) = 0.5.*(RHO(i,:) + ones(size(RHO(i,:)))); % assume density is 1 at farfield
        else
            RHO_Tavg(i,:) = 0.5.*(RHO(i,:) + RHO(i+1,:));
        end
    end
    
    RHO_Ravg = zeros(size(RR));
    for i = 1:size(RHO_Ravg,2) % loop through Thetas in original grid
        if i == 1 % Apply Neumann condition
           i1 = i;
           i_1 = i1;
        elseif i == size(RHO_Ravg,2) % Apply Neumann Condition
            i1 = i-1;
            i_1 = i1;
        else
            i1 = i;
            i_1 = i-1;
        end
        RHO_Ravg(:,i) = 0.5.*(RHO(:,i1) + RHO(:,i_1)); 
    end
    
    for i = 1:length(T_vals) % loop through theta        
        for j = 1:(length(r_vals)-1) % loop through radius
            % 5 different control indexes for i and j directions            
            i0 = i;
            j0 = j;
            % Boundary Conditions for Theta
            if i == 1
                i1 = i+1;
                i_1 = i1;
            elseif i == length(T_vals)
                i1 = i-1;
                i_1 = i-1;
            else
                i1 = i+1;
                i_1 = i-1;
            end
            
            % Boundary Conditions for Density
            if i == 1 % Apply Neumann condition
                i1r = i;
                i_1r = i1r;
            elseif i == size(RHO_Ravg,2) % Apply Neumann Condition
                i1r = i-1;
                i_1r = i1r;
            else
                i1r = i;
                i_1r = i-1;
            end

            % Boundary Conditions for Radius
            if j == 1
                j1 = j+1;
                j_1 = j1;
            else
                j1 = j+1;
                j_1 = j-1;
                RR_j0 = 0.5*(RR(j_1, i0) + RR(j0,i0));
            end
            
            % Radius Averages
            RR_j1 = 0.5*(RR(j1, i0) + RR(j0,i0));
            
            if j == 1
                PHI_RR = (RR_j1*RHO_Ravg(j1,i0)*(PHI(j1,i0,iter) - PHI(j0,i0,iter))/dr)/(RR(j0,i0)*dr);
            else
                PHI_RR = (RR_j1*RHO_Ravg(j1,i0)*(PHI(j1,i0,iter) - PHI(j0,i0,iter))/dr - RR_j0*RHO_Ravg(j0,i0)*(PHI(j0,i0,iter) - PHI(j_1,i0,iter))/dr)/(RR(j0,i0)*dr);
            end
            
            PHI_TT = (1/RR(j0,i0)^2)*(1/dT)*(RHO_Tavg(j0, i1r)*(PHI(j0,i1,iter) - PHI(j0,i0,iter))/dT - RHO_Tavg(j0, i_1r)*(PHI(j0,i0,iter) - PHI(j0,i_1,iter))/dT);
            
            PHI(j0,i0,iter+1) = ((0.5*alpha/dt - 1/(dt^2))*PHI(j0,i0,iter-1) + 2/(dt^2)*PHI(j0,i0,iter) + PHI_RR + PHI_TT)/(1/dt^2 + 0.5*alpha/dt);
            
        end
    end
    
    PHI(end,:,iter+1) = BC.dirichlet.PHI_II; % apply BC at outer boundary
    
    % Run Residual Check
    difference = abs(PHI(:,:,iter+1) - PHI(:,:,iter));
    [res(iter), ind(iter)] = max(difference(:));
end

%% Post Process

figure(1);semilogy(1:iter, res);

% plot density
figure();contourf(XX,YY,RHO_Ravg, 50)

figure(); % cp plots
contourf(XX, YY, 1-(q2_ij), 50); %./((RR.*cos(TT)).^2)
title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal

figure();
plot(TT(1,:)*180/pi, 1 - q2_ij(1,:));
xlabel('\theta');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title('C_p on surface of Cylinder');
set(gca, 'Xdir', 'reverse');

figure(); % field potential
contourf(XX, YY, PHI(:,:,end), 50);
title('Field Potential');
colorbar('eastoutside');
axis equal

figure(); % x-dir velocity plots
contourf(XX, YY, U_n, 50); %./((RR.*cos(TT)).^2)
title('PHI_\theta velocity');
colorbar('eastoutside');

figure();
contourf(XX, YY, sqrt(M2_ij), 50);
title('Mach Number');
colorbar('eastoutside');