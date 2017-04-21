clc;
clear;
close all;

%% SIM CONTROL PARAMS - GRID INITIALIZATION

% Step Sizes
dx = 0.1;
dy = dx;
dt = 0.1*dx;

alpha = 10.0;

% Airfoil Dimensions
tau = 0.05;
chord = 1.0;

% Field Axis Values
y_max = 5;
x_max = 6.0;
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
M0 = 0.0;
visc_on = 0;

%% SIM CONTROL VARIABLE INITIALIZATION

PHI = zeros([size(XX), 1]);
PHI(:,:,1) = XX;
% PHI(:,:,2) = PHI(:,:,1);
[n_r, n_T] = size(XX);
res = 1;
ind = 0;
iter = 0.0;
tol = 1e-4;

% Boundary Conditions
BC.Vy_II = zeros(size(YY(end,:)));
BC.Vx_II = ones(size(XX(end,:)));
BC.Vx_I = ones(size(XX(:,1)));
% BC.PHI_II = XX(end,:);
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
    PHI(:,:,iter+1) = PHI(:,:,iter); % setup iter+1th iteration step
    
    for ii = 2:size(PHI,2) % loop through x-values
        % Calculate Densities in front and behind
        u1 = (PHI(:,ii+1) - PHI(:,ii)) / dx; % next velocity
        u_1 = (PHI(:,ii) - PHI(:,ii-1)) / dx; % previous velocity
        
        % Initial Density Calcs
        rho1 = (1 - 0.5*(gam-1)*M0^2.*(u1.^2 -1)).^(1/(gam-1));
        rho_1 = (1 - 0.5*(gam-1)*M0^2.*(u_1.^2 -1)).^(1/(gam-1));
        
        
        
    end

end

%% Post Process

% figure();contourf(XX,YY,RHO_Ravg, 50)

figure(); % cp plots
contourf(XX, YY, 1-(U_n.^2), 50); %./((RR.*cos(TT)).^2)
title('Pressure Coefficient Contours');
colorbar('eastoutside');
axis equal

figure();
plot(XX(1,:), 1 - U_n(1,:).^2);
xlabel('\theta');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title('C_p on surface of Airfoil');
set(gca, 'Xdir', 'reverse');

figure(); % field potential
contourf(XX, YY, PHI(:,:,end), 50);
title('Field Potential');
colorbar('eastoutside');
axis equal

figure(); % x-dir velocity plots
contourf(XX, YY, U_n, 50); %./((RR.*cos(TT)).^2)
title('\Phi_x velocity');
colorbar('eastoutside');

figure();
contourf(XX, YY, M_ij, 50);
title('Mach Number');
colorbar('eastoutside');