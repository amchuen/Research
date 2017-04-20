clc;
clear;
close all;

%% SIM CONTROL PARAMS - GRID INITIALIZATION

% Step Sizes
dx = 1/40;
dy = dx;
dt = 0.1*dx;

alpha = 5.0;

% Airfoil Dimensions
tau = 0.05;
chord = 1.0;

% Field Axis Values
y_max = 10;
x_max = 5.0;
x_vals = (-4*dx):dx:x_max;
y_vals = 0:dy:y_max;
[XX, YY] = meshgrid(x_vals, y_vals);

% Body Values
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

% Fluid Params
gam = 1.4; % heat 
M0 = 0.25;
visc = 0.0;

%% SIM CONTROL VARIABLE INITIALIZATION

PHI = zeros([size(XX), 2]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
[n_r, n_T] = size(XX);
res = 1;
ind = 0;
iter = 1.0;
tol = 1e-4;

% Boundary Conditions
BC.Vy_II = zeros(size(YY(end,:)));
BC.Vx_II = ones(size(XX(end,:)));
BC.Vx_I = ones(size(XX(:,1)));
BC.PHI_II = XX(end,:);
BC.PHI_I = XX(:,1);

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
    
    % Apply boundary conditions
    PHI(end,:,iter+1) = BC.PHI_II;% farfield condition
    PHI(:,1,iter+1) = BC.PHI_I; % inlet condition
    PHI(:,2,iter+1) = XX(:,2); % inlet condition

%     % Calculate Density from previous sweep through field (lagging?)
%     U_n = ones(size(XX)); % PHI_X... only need velocity in X-dir
%     
%     for i = 2:(size(XX, 2)-1) % loop through x-vals
%         U_n(:,i) = (PHI(:,i+1,iter) - PHI(:,i-1,iter))./(2.*dx);
%         % determine velocity at the far-field end?
%     end
%     U_n(:,1) = BC.Vx_I;
    
    % Initialize Density calculations
    rho_ij = ones([size(XX,1), size(XX,2)-1]);
    u_avg = rho_ij;
    
    for i = 2:(n_r) % looping through y-dir points
        for ii = 1:(size(u_avg,2)) % looping through x-dir points
                % Will have 1 less point than grid b/c of offset
                % Neumann enforced for density by setting end grid-points to be
                % same if reflecte across boundary

                % Calculate average velocities... take average of four
                % surrounding points of elements
    %             u_avg(i,ii) = mean([U_n(i, ii), U_n(i-1, ii), U_n(i-1, ii+1), U_n(i, ii+1)]);
            u_avg(:,ii) = 0.5*((PHI(i,ii+1,iter) - PHI(i,ii,iter))/dx + (PHI(i-1,ii+1,iter) - PHI(i-1,ii,iter))/dx);
        end
    end
    
    % Apply Neumann condition on velocities
    u_avg(1,:) = u_avg(2,:);
    
    a2_avg = (1./M0.^2)-0.5*(gam-1).*(u_avg.^2 - 1);
    eps_ij = max(zeros(size(rho_ij)), 1 - 0.9.*(a2_avg./(u_avg.^2)));
    
    if any(any((1 - 0.5*(gam-1).*M0.^2.*(u_avg.^2 - 1))<0))
       fprintf('Non-real result for density!\n'); 
    end
    
    % Calculate uncorrected density 
    rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(u_avg.^2 - 1)).^(1./(gam-1));
    
%     if rho_ij(1,:) ~= rho_ij(2,:) % check if Neumann condition satisfied by cylinder surface
%         fprintf('Neumann condition not satisfied for density! Please check density before continuing.\n');
%     end
    
    % Calculate viscosity corrections
%     rho_ds = [rho_ij(:,1)-ones(size(rho_ij(:,1))), diff(rho_ij,1,2)];
    rho_ds = zeros(size(rho_ij));
    for i = 1:(size(rho_ij,2)) % loop through x-dir
        if i ==1 % Dirichlet BC
            rho_ds(:,i) = rho_ij(:,i)-ones(size(rho_ij(:,i)));
        else
            rho_ds(:,i) = rho_ij(:,i)-rho_ij(:,i-1);
        end
    end
    
    RHO = rho_ij - eps_ij.*rho_ds;
    
    RHO_Xavg = RHO;%zeros(size(RHO));
%     for i = 1:size(RHO,1) % loop through y-dir
%         if i == size(RHO,1)
%             RHO_Xavg(i,:) = 0.5.*(RHO(i,:) + ones(size(RHO(i,:)))); % assume density is 1 at farfield
%         else
%             RHO_Xavg(i,:) = 0.5.*(RHO(i,:) + RHO(i+1,:)); % assumes density grid is shifted half down from regular grid
%         end
%     end
    
    for i = 2:(size(PHI,2)-1) % march along x-dir to calculate one step ahead, do not need to calculate final wall
        for j = 1:(size(PHI,1)-1) % loop through y_dir, do not need to calculate for top BC
            if j == 1
%                 % Apply Neumann Condition... either air or body
%                 if (XX(j,i) >= 0) && (XX(j,i) <=1) % body condition
%                     PHI_Y_1 = YY
%                 else
%                     PHI_Y_1 = 0; % air condition
%                 end
                PHI_Y_1 = dyBdx(i);
            else % everywhere else not on body
                PHI_Y_1 = (PHI(j,i,iter+1) - PHI(j-1,i,iter+1))/dy;
            end
            RHO_i1 = RHO_Xavg(j,i);
            if i == 2
                RHO_i_1 = 1; % assumes free-stream density of 1
            else
                RHO_i_1 = RHO_Xavg(j,i-1);
            end
            
            PHI_YY = ((PHI(j+1,i,iter+1)-PHI(j,i,iter+1))/dy - PHI_Y_1)/dy;
            PHI(j,i+1,iter+1) = PHI(j,i,iter+1) - (dx^2 / RHO(j,i)).*(PHI_YY) + (RHO_i_1/RHO_i1)*(PHI(j,i,iter+1) - PHI(j,i-1,iter+1));
        end
        
        if XX(1,i) >= 0
            close all;
            figure();plot(y_vals, PHI(:,i,end))
            fprintf('Please check calculations here!\n'); 
        end
    end
    
    % Run Residual Check
    difference = abs(PHI(:,:,iter+1) - PHI(:,:,iter));
    [res(iter), ind(iter)] = max(difference(:));
end