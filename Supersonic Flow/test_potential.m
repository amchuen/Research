% clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% % Field Axis Values
% y_max = 20;
% x_max = 10 + 20*GR.dx;
% x_min = -10-39*GR.dx; %(-19*dx);
% GR.x_vals = x_min:GR.dx:x_max;
% GR.y_vals = 0:GR.dy:y_max;
% [GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

% Field Axis Values
y_max = 7;
x_max = 2 + 20*GR.dx;
x_min = -2-39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
FL.M0 = 1.401;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 1.3;
CT.tol = 1e-5;
CT.alpha = 50;
CT.dt = 0.05;
CT.mom_corr = 0;
CT.t_rho = 1;
CT.enforce_phi = 1; % enforces boundary condition for potential
CT.enforce_vel = 1;
CT.enforce_CFL = 1;
CT.phi_xx = 1; % use zeros at exit BC

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.05;
x_vals = GR.x_vals;
dx = GR.dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

BC.Vy_II = zeros(size(GR.YY(end,:)));
BC.Vx_II = ones(size(GR.XX(end,:)));
BC.Vx_I = ones(size(GR.XX(:,1)));
BC.PHI_II = GR.XX(end,:);
BC.PHI_I = GR.XX(:,1);
BC.dyBdx = dyBdx;

%% RUN CALCULATION

close all;

PHI = zeros([size(GR.XX), 2]);
PHI(:,:,1) = GR.XX;
PHI(:,:,2) = PHI(:,:,1);
PHI_new = PHI(:,:,2);
rho_ij = ones(size(GR.XX(:,1:(end-1))));
rho_ds = zeros(size(rho_ij));
res = 1;
visc_on = 0;
ind = [];

%% START SOLUTION

while (res(end) > CT.tol)|| (length(res) < 300) % iterate through time

    % Reorganize three-level scheme
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI_new;

    % Apply boundary conditions
    PHI_new(:,1) = BC.PHI_I; % inlet condition
    
    % Velocity -> first order
    phiX_w = [BC.Vx_I, (PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./GR.dx];
    if FL.M0>1
        phiX_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./GR.dx, (PHI(:,end, 2) - PHI(:,end-1,2))./GR.dx];
    else
        phiX_e = [(PHI(:, 2:end, 2) - PHI(:,1:(end-1),2))./GR.dx, ones(size(PHI(:,end,2)))];
    end
    phiY_n = [(diff(PHI(:,:,2),1,1)./GR.dy); (GR.XX(end,:) - PHI(end,:,2))./GR.dy];
    phiY_s = [BC.dyBdx; (diff(PHI(:,:,2),1,1)./GR.dy)]; % enforce irrotationality at the body

    %% Calculate local viscosity at the nodes
    phiX = 0.5.*(phiX_w + phiX_e);
    phiY = 0.5.*(phiY_n + phiY_s);
    if max(abs(phiX(:,1) - BC.Vx_I)) > CT.tol && CT.enforce_Vel
        phiX(:,1) = BC.Vx_I;
    end
    
    if FL.M0 > 1
       phiX(:,end) =  phiX(:,end-1);
    end

    a2_ij = (1./FL.M0.^2)-0.5*(FL.gam-1).*(phiX.^2 - 1);
    eps_ij = max(zeros(size(phiX)), 1 - 1.0.*(a2_ij./(phiX.^2)));

    if any(any(eps_ij ~= 0)) && (visc_on ==0) % viscosity not being used and then turned on
        visc_on = 1;
        fprintf('Viscosity model activated! Iteration: %i\n', length(res));        
    elseif (visc_on == 1) && all(all(eps_ij == 0))
        visc_on = 0;
        fprintf('Viscosity model turned off! Iteration: %i\n', length(res));
    end
    
    % Check CFL condition
    CFL_i = (max(abs(phiX(:)))./GR.dx + max(abs(phiY_s(:)))./GR.dy)*CT.dt;
    
    if CFL_i >= 1.0 && CT.enforce_CFL
       fprintf('CFL condition not met @ timestep %i! Decreasing time steps!\n', length(res));
       CT.dt = CT.dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', CT.dt);
    end

    %% Calculate Density using between grid scheme in X-dir
    phi_t = (PHI(:,:,2) - PHI(:,:,1))./CT.dt;
    phi_t_avg = 0.5.*(phi_t(:,1:(end-1)) + phi_t(:,2:end)); % average phi_t for three level scheme
    rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(phiX_e(:,1:(end-1)).^2 + 2.*phi_t_avg - 1)).^(1./(FL.gam-1)); % do I need the time correction?

    if any(imag(rho_ij(:)) ~= 0)
%         dbstop;
        test = imag(rho_ij);
        [y_ind, x_ind] = find(test ~= 0);
        fprintf('Non-real result for density!\n'); 
        fprintf('Check Nodes...\n');
        fprintf('Row %i, Col %i \n', y_ind, x_ind);
        fprintf('Limitting Density to Zero! \n');
        rho_ij(y_ind, x_ind) = 0;
    end

    % Calculate viscosity corrections
    rho_ds(:,1) = 0.5*(rho_ij(:,2)-ones(size(rho_ij(:,1)))) - 0.5.*sign(phiX_e(:,1)).*(rho_ij(:,2)-2.*rho_ij(:,1) + ones(size(rho_ij(:,1))));
    rho_ds(:,end) = 0.5*(-rho_ij(:,end-1)+ones(size(rho_ij(:,end)))) - 0.5.*sign(phiX_e(:,end-1)).*(rho_ij(:,end-1)-2.*rho_ij(:,end) + ones(size(rho_ij(:,end))));
    rho_ds(:,2:(end-1)) = 0.5*(rho_ij(:,3:end)-rho_ij(:,1:(end-2))) - 0.5.*sign(phiX_e(:,2:(end-2))).*(rho_ij(:,3:end)-2.*rho_ij(:,2:(end-1)) + rho_ij(:,1:(end-2)));

    if (FL.M0 == 0) && any(any(rho_ij ~= 1))
        fprintf('Incompressible Case not satisfied!\n');
    end
    
    %% Three-Level Scheme
%     if FL.M0 < 1
        % Calculate Density in X-Dir
        RHO_W = rho_ij(:,1:(end-1)) - CT.v_coeff.*eps_ij(:,2:(end-1)).*rho_ds(:,1:(end-1));
        RHO_E = rho_ij(:,2:end) - CT.v_coeff.*eps_ij(:,2:(end-1)).*rho_ds(:,2:end); % use local viscosity at current node

        % Calculate Laplacian coupled with density
        PHI_XX = (RHO_E.*phiX_e(:,2:(end-1)) - RHO_W.*phiX_w(:,2:(end-1)))./GR.dx;
        PHI_YY = (phiY_n(:,2:(end-1)) - phiY_s(:,2:(end-1)))./GR.dy;

        % Include Momentum Correction
        RHO = 0.5*(RHO_E + RHO_W);
        if CT.mom_corr
            CORR = ((RHO.^(FL.gam - 1) - 1)./((FL.gam-1)*FL.M0^2) + 0.5.*(phiX(:,2:(end-1)).^2 - 1));
        else
            CORR = zeros(size(RHO));
        end

        % Extract Solution
        PHI_new(:,2:(end-1)) = (PHI_XX + PHI_YY - CT.alpha.*CORR + 2/(CT.dt^2)*PHI(:,2:(end-1),2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,2:(end-1),1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);

    if FL.M0 > 1 % supersonic, use Hyperbolic Solver?
        
        if CT.phi_xx == 1 % Extrapolate from upstream results
            RHO_i1 = rho_ij(:,end) - CT.v_coeff.*eps_ij(:,end).*rho_ds(:,end); % use local viscosity at current node
            RHO_i_1 = rho_ij(:,end-1) - CT.v_coeff.*eps_ij(:,end).*rho_ds(:,end-1);
            PHI_XX = (RHO_i1.*phiX_e(:,end-1) - RHO_i_1.*phiX_w(:,end-1))./GR.dx;
        elseif CT.phi_xx == 2 % Explicit assumption of PHI_XX
            PHI_XX = zeros(size(GR.XX(:,end)));
        end

        PHI_YY = (phiY_n(:,end) - phiY_s(:,end))./GR.dy;

        if CT.mom_corr && (CT.phi_xx == 1)
            RHO_end = 0.5*(RHO_i1 + RHO_i_1);
            CORR_end = ((RHO_end.^(FL.gam - 1) - 1)./((FL.gam-1)*FL.M0^2) + 0.5.*(phiX(:,end).^2 - 1));
        else
            CORR_end = zeros(size(PHI_YY));
        end

        PHI_new(:,end) = (PHI_XX + PHI_YY - CT.alpha.*CORR_end + 2/(CT.dt^2)*PHI(:,end,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,end,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
    end
%     elseif FL.M0 > 1 % supersonic, use Hyperbolic Solver?
% 
%         for i = 3:(size(PHI,2)) % march along x-dir to calculate one step ahead, do not need to calculate initial inlet
%             if i == size(PHI,2)
%                 RHO_i1 = rho_ij(:,i-1) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-1); % use local viscosity at current node
%                 RHO_i_1 = rho_ij(:,i-2) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-2);
% 
%                 PHI_XX = (RHO_i1.*(PHI(:,i,2) - PHI(:,i-1,2))./GR.dx - RHO_i_1.*(PHI(:,i-1,2) - PHI(:,i-2,2))./GR.dx)./GR.dx;
%                 
% %                 
% %                 PHI_YY = sys_yy*PHI(:,i,2) - (1/GR.dy).*[BC.dyBdx(i); zeros(length(GR.y_vals)-1, 1)];
%                 PHI_YY = [((PHI(2:end,i,2) - PHI(1:end-1,i,2))./GR.dy - [BC.dyBdx(i); (PHI(2:end-1,i,2) - PHI(1:end-2,i,2))./GR.dy])./GR.dy; 0];
% 
%                 PHI_new(:,i) = (PHI_XX + PHI_YY + 2/(CT.dt^2)*PHI(:,i,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,i,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
% %                 PHI_new(:,i) = PHI_new(:,i-1);
%             else % start from line three
% 
%                 RHO_i1 = rho_ij(:,i) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i); % use local viscosity at current node
%                 RHO_i_1 = rho_ij(:,i-1) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-1);
% 
%                 PHI_XX = (RHO_i1.*(PHI(:,i+1,2) - PHI(:,i,2))./GR.dx - RHO_i_1.*(PHI(:,i,2) - PHI(:,i-1,2))./GR.dx)./GR.dx;
% 
% %                 PHI_YY = sys_yy*PHI(:,i,2) - (1/GR.dy).*[BC.dyBdx(i); zeros(length(GR.y_vals)-1, 1)];
%                 PHI_YY = [((PHI(2:end,i,2) - PHI(1:end-1,i,2))./GR.dy - [BC.dyBdx(i); (PHI(2:end-1,i,2) - PHI(1:end-2,i,2))./GR.dy])./GR.dy; 0];
% 
%                 PHI_new(:,i) = (PHI_XX + PHI_YY + 2/(CT.dt^2)*PHI(:,i,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,i,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
%             end
%         end
%         PHI_new(end,:) = BC.PHI_II;
%     end
    
    % Check if upper far-field BC's are met
    if (max(abs(PHI_new(end,:) - BC.PHI_II)) > CT.tol) && CT.enforce_phi
        PHI_new(end,:) = BC.PHI_II;
    end

    % Run Residual Check
    difference = abs(PHI_new(:,:) - PHI(:,:,2));
    [res(end+1), ind(end+1)] = max(difference(:));

    % Initialize next time step
    if (length(res) > 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
        figure(1);semilogy(1:length(res), res);
    end

end

%% POST PROCESS

phiX = ones(size(GR.XX)); % PHI_X... only need velocity in X-dir
PHI(:,:,2) = PHI_new;

for i = 2:(size(GR.XX, 2)-1) % loop through x-vals
    phiX(:,i) = (PHI(:,i+1,end) - PHI(:,i-1,end))./(2.*GR.dx);
    % determine velocity at the far-field end?
end
phiX(:,1) = BC.Vx_I;

if FL.M0 > 1
   phiX(:,end) =  phiX(:,end-1);
end

% Calculate uncorrected density 
%     rho_ij = (1 - 0.5*(FL.gam-1).*M  0.^2.*(u_avg.^2 - 1)).^(1./(FL.gam-1));
rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(phiX.^2 - 1)).^(1./(FL.gam-1));
a2_avg = (1./FL.M0.^2)-0.5*(FL.gam-1).*(phiX.^2 - 1);
M_ij = phiX./sqrt(a2_avg);

% Assign to output
OUT.PHI = PHI;
OUT.U_n = phiX;
OUT.a_l = sqrt(a2_avg);
OUT.rho_ij = rho_ij;
OUT.M_ij = M_ij;
OUT.res = res;