clc;
clear;
close all;

%% SIM CONTROL PARAMS - GRID INITIALIZATION

% Step Sizes
dx = 0.1;
dy = dx;
dt = 0.1*dx;

alpha = 25;

% Airfoil Dimensions
tau = 0.05;
chord = 1.0;

% Field Axis Values
y_max = 5;
x_max = 6.0;
x_vals = (-4):dx:x_max;
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
M0 = 1.1;
visc_on = 0;

%% SIM CONTROL VARIABLE INITIALIZATION

PHI = zeros([size(XX), 2]);
PHI(:,:,1) = XX;
PHI(:,:,2) = PHI(:,:,1);
[n_r, n_T] = size(XX);
res = 1;
ind = 0;
iter = 1.0;
tol = 5e-5;

% Boundary Conditions
BC.Vy_II = zeros(size(YY(end,:)));
BC.Vx_II = ones(size(XX(end,:)));
BC.Vx_I = ones(size(XX(:,1)));
% BC.PHI_II = XX(end,:);
BC.PHI_I = XX(:,1);

%% START SOLUTION

while (res(end) > tol) % iterate through time
    
    % Initialize next time step
    if (iter > 500) && (mod(iter, 500) == 0)
        fprintf('Iteration Ct: %i\n', iter);
        fprintf('Current Residual: %0.5f\n', res(end));
        figure(1);semilogy(1:iter, res);
    end
    iter = iter + 1; % use this to call nth time step
    PHI(:,:,iter+1) = PHI(:,:,iter); % setup n+1th time step 
    
    % Apply boundary conditions
    PHI(:,1,iter+1) = BC.PHI_I; % inlet condition

    % Initialize Density calculations
    rho_ij = ones([size(XX,1), size(XX,2)-1]);
    u_avg = diff(PHI(:,:,iter)')'./dx;
    
    % Calculate Local Viscosity
    a2_avg = (1./M0.^2)-0.5*(gam-1).*(u_avg.^2 - 1);
    eps_ij = max(zeros(size(rho_ij)), 1 - 0.9.*(a2_avg./(u_avg.^2)));
    
    % Check if viscosity is being used
    if any(any(eps_ij ~= 0)) && (visc_on ==0) % viscosity not being used and then turned on
        visc_on = 1;
        fprintf('Viscosity model activated! Iteration: %i\n', iter);        
    elseif (visc_on == 1) && all(all(eps_ij == 0))
        visc_on = 0;
        fprintf('Viscosity model turned off! Iteration: %i\n', iter);
    end
    
    rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(u_avg.^2 - 1)).^(1./(gam-1));
    
    if any(any((1 - 0.5*(gam-1).*M0.^2.*(u_avg.^2 - 1))<0))
       fprintf('Non-real result for density!\n'); 
    end
    
    % Calculate viscosity corrections
    rho_ds = zeros(size(rho_ij));
    for i = 1:(size(rho_ij,2)) % loop through x-dir
        if i ==1 % Dirichlet BC
            rho_ds(:,i) = rho_ij(:,i)-ones(size(rho_ij(:,i)));
        else
            rho_ds(:,i) = rho_ij(:,i)-rho_ij(:,i-1);
        end
    end
    
    RHO = rho_ij - eps_ij.*rho_ds;
        
    if (M0 == 0) && any(any(RHO ~= 1))
        fprintf('Incompressible Case not satisfied!\n');
    end
    
    for i = 2:(size(PHI,2)-1) % march along x-dir to calculate one step ahead, do not need to calculate initial inlet
        for j = 1:(size(PHI,1)-1) % loop through y_dir, do not need to calculate for top BC
            if j == 1 % applies either body or foil
                PHI_Y_1 = dyBdx(i);
            else % everywhere else not on body
                PHI_Y_1 = (PHI(j,i,iter+1) - PHI(j-1,i,iter+1))/dy;
            end
            
            if i == size(PHI,2)
                RHO_i1 = 1;
            else
                RHO_i1 = RHO(j,i);
            end
            RHO_i_1 = RHO(j,i-1);
            
            PHI_YY = ((PHI(j+1,i,iter)-PHI(j,i,iter))/dy - PHI_Y_1)/dy;
            PHI_XX = (RHO_i1*(PHI(j,i+1,iter) - PHI(j,i,iter))/dx - RHO_i_1*(PHI(j,i,iter) - PHI(j,i-1,iter))/dx)/dx;
            PHI(j,i,iter+1) = (PHI_XX + PHI_YY + 2/(dt^2)*PHI(j,i,iter) - (1/dt^2 - 0.5*alpha/dt)*PHI(j,i,iter-1))/(1/dt^2 + 0.5*alpha/dt);
        end
    end
    
    % Run Residual Check
    difference = abs(PHI(:,:,iter+1) - PHI(:,:,iter));
    [res(iter), ind(iter)] = max(difference(:));
    
end

%% Post Process

% Calculate Density from previous sweep through field (lagging?)
U_n = ones(size(XX)); % PHI_X... only need velocity in X-dir

for i = 2:(size(XX, 2)-1) % loop through x-vals
    U_n(:,i) = (PHI(:,i+1,iter) - PHI(:,i-1,iter))./(2.*dx);
    % determine velocity at the far-field end?
end
U_n(:,1) = BC.Vx_I;

% Calculate uncorrected density 
%     rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(u_avg.^2 - 1)).^(1./(gam-1));
rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(U_n.^2 - 1)).^(1./(gam-1));
a2_avg = (1./M0.^2)-0.5*(gam-1).*(U_n.^2 - 1);
M_ij = U_n./sqrt(a2_avg);


folderName = ['M_' num2str(M0)];

if ~exist([pwd '\airfoil\' folderName], 'dir')
    mkdir([pwd '\airfoil\' folderName]);
end

close all;
figure(1);semilogy(1:iter, res);
title('Residual Plot');
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd '\airfoil\' folderName '\residual_plot.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,rho_ij, 50)
title('Density (Normalized)');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\airfoil\' folderName '\density.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-(U_n.^2), 50); %./((RR.*cos(TT)).^2)
title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\airfoil\' folderName '\cp_contour.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\cp_contour']);

figure();
plot(XX(1,:), 1 - U_n(1,:).^2);
xlabel('\theta');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title('C_p on surface of airfoil');
set(gca, 'Ydir', 'reverse');
saveas(gcf, [pwd '\airfoil\' folderName '\cp_surf.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\cp_surf']);

figure(); % field potential
contourf(XX, YY, PHI(:,:,end), 50);
title('Field Potential, \Phi');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\airfoil\' folderName '\phi_pot.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, U_n, 50); %./((RR.*cos(TT)).^2)
title('\Phi_\theta velocity');
colorbar('eastoutside');
saveas(gcf, [pwd '\airfoil\' folderName '\phi_theta.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\phi_theta']);

figure();
contourf(XX, YY, M_ij, 50);
title('Mach Number');
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\airfoil\' folderName '\mach.png']);
saveas(gcf, [pwd '\airfoil\' folderName '\mach']);