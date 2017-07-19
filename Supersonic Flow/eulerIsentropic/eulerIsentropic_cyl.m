clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

case_name = 'cylinder';

% Define Grid
dT = 0.01*pi;
dr = 0.05;

% Field Axis Values
r_max = 25;
r_cyl = 0.5; 
x_max = 30; %10 + 20*dr;
x_min = 0.0; %-39*dr; %(-19*dr);
T_vals = 0.0:dT:pi;%x_min:dT:x_max;
R_vals = r_cyl:dr:r_max;
[TT, RR] = meshgrid(T_vals, R_vals);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.25;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.000075; % spatial diffusion term
eps_t = 0.000075; % time diffusion term
tol = 1e-3;
dt = 0.01^2;
iter_min = 300;
CFL_on = 1;

%% INITIALIZE

% Three level scheme
RHO.fv = ones([size(TT), 3]); % density
AA.fv = repmat(-(1 + (r_cyl^2)./(RR.^2)).*sin(TT), 1, 1, 3); % rho * u (theta)
BB.fv = repmat(cos(TT).*(1 - (r_cyl^2)./((RR).^2)),1,1,3); % rho * v (radial vel)
PP.fv = RHO.fv(:,:,2).^(gam) ./ (gam .* M0.^2); % pressure

%% BodT Values - Cylinder
tau = 2;%0.1;
m_x = tand(8); % dT/dr
YY_B = [zeros(size(T_vals(T_vals <6))), ...
        tau.*ones(size(T_vals(T_vals >= 6)))];
dTBdr = zeros(size(YY_B));


%% Boundary Conditions

% Inlet - West (Symmetry)
RHO.BC.Wx = zeros(size(TT(:,end))); % neumann condition
AA.BC.Wx = RHO.BC.Wx; % neumann condition
BB.BC.Wx = RHO.BC.Wx; % neumann condition
AB_R.BC.Wx = RHO.BC.Wx; % Neumann condition
PP.BC.Wx = RHO.BC.Wx; % Neumann condition

% Outlet - East (Symmetry)
RHO.BC.Ex = zeros(size(TT(:,end))); % neumann condition
AA.BC.Ex = RHO.BC.Ex; % neumann condition
BB.BC.Ex = RHO.BC.Ex; % neumann condition
AB_R.BC.Ex = RHO.BC.Ex; % Neumann condition
PP.BC.Ex = RHO.BC.Ex; % Neumann condition
%(RHO.fv(:,end,2).^(gam))./(gam .* M0.^2); % dirichlet BC, update with time

% Far-field - North
RHO.BC.N = ones(size(TT(end,:))); % dirichlet BC
BB.BC.N = -(1 + (r_cyl^2)./((RR(end,:)).^2)).*sin(TT(end,:)); % dirichlet BC
AA.BC.N = cos(TT(end,:)).*(1 - (r_cyl^2)./((RR(end,:)).^2)); % dirichlet BC
PP.BC.N = (RHO.BC.N.^(gam))./(gam .* M0.^2); % dirichlet BC
AB_R.BC.N = AA.BC.N .* BB.BC.N ./ RHO.BC.N; % dirichlet BC

% Wall/Symmetry - South
RHO.BC.Sy = zeros(size(TT(1,:))); % neumann condition
BB.BC.Sy = zeros(size(TT(1,:))); % neumann condition
AA.BC.S = RHO.fv(1,:,2).*dTBdr; % dirichlet condition... update with time
PP.BC.Sy = zeros(size(TT(1,:))); % neumann condition
AB_R.BC.S = BB.BC.S .* AA.fv(1,:,2) ./ RHO.fv(1,:,2); % dirichlet condition... update with time

%% START SOLUTION

res = [1, 1, 1];
tic;
while (max(res(end, :)) > tol*max(res(:)))|| (size(res,1) < iter_min) % iterate through time
    
    % Update Field Values
    RHO.fv(:,:,1) = RHO.fv(:,:,2);
    RHO.fv(:,:,2) = RHO.fv(:,:,3);
    AA.fv(:,:,1) = AA.fv(:,:,2);
    AA.fv(:,:,2) = AA.fv(:,:,3);
    BB.fv(:,:,1) = BB.fv(:,:,2);
    BB.fv(:,:,2) = BB.fv(:,:,3);
    PP.fv = RHO.fv(:,:,2).^(gam) ./ (gam .* M0.^2); % pressure
    
    % Update Boundary Conditions
%     BB.BC.S = RHO.fv(1,:,2).*dTBdr; % dirichlet condition... update with time
    AB_R.S = BB.BC.S .* AA.fv(1,:,2) ./ RHO.fv(1,:,2); % dirichlet condition... update with time
    
    % Check CFL conditions
    Ux = (AA.fv(:,:,2)./RHO.fv(:,:,2));
    Vy = (BB.fv(:,:,2)./RHO.fv(:,:,2));
    CFL_i = (max(abs(Ux(:)))./dT + max(abs(Vy(:)))./dr)*dt;

    if CFL_i >= 1.0 && CFL_on
       fprintf('CFL condition not met!\n');
%            if CFL_on
       fprintf('Decreasing time steps!\n');
       dt = dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
%            end
%     elseif CFL_i < 0.8
%         dt = dt / CFL_i;
    end
    
    % Calculate Derivatives
    [Ar, ~, ~] = grad_f(RR.*AA.fv(:,:,2), 1, dr, AA.BC, 1);
    [BT, ~, ~] = grad_f(BB.fv(:,:,2), 2, RR.*dT, BB.BC, 1);
    [AT, ~, ~] = grad_f(AA.fv(:,:,2), 2, RR.*dT, AA.BC, 1);
    [RHOT, ~, ~] = grad_f(RHO.fv(:,:,2), 2, RR.*dT, RHO.BC, 1);
    [Br, ~, ~] = grad_f(BB.fv(:,:,2), 1, dr, BB.BC, 1);
    [RHOr, ~, ~] = grad_f(RHO.fv(:,:,2), 1, dr, RHO.BC, 1);
    
    % Calculate new RHO
	laplace_RHO = laplace_f(RHO.fv(:,:,2), dr, dT, RHO.BC, 1); % note -> need to figure out how to modify radial term for polar coordinates
	
    RHO.fv(:,:,3) = (eps_s.*laplace_RHO - Ar./RR - BT - (eps_t./(dt^2) - 0.5./dt).*RHO.fv(:,:,1) + 2.*eps_t.*RHO.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    % Enforce BC's
    RHO.fv(:,1,3) = RHO.BC.W;
    RHO.fv(end,:,3) = RHO.BC.N;    
    
    % Calculate new A
%     [A2_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).^2)./RHO.fv(:,:,2), 2, dr, AA.BC, 1);
%     [AB_rho_y, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 1, dT, AB_R.BC, 1);
    A2_rho_x = 2.*AA.fv(:,:,2).*Ar./RHO.fv(:,:,2) - (AA.fv(:,:,2).^2 .* RHOr)./(RHO.fv(:,:,2).^2);
    AB_rho_y = (AA.fv(:,:,2).*BT./RHO.fv(:,:,2)) + (BB.fv(:,:,2).*AT)./(RHO.fv(:,:,2)) - (AA.fv(:,:,2).*BB.fv(:,:,2).*RHOT)./(RHO.fv(:,:,2).^2);
    laplace_A = laplace_f(AA.fv(:,:,2), dr, dT, AA.BC, 1);
    [P_x, ~, ~] = grad_f(PP.fv, 2, dr, PP.BC, 1);
    
    AA.fv(:,:,3) = (eps_s.*laplace_A - P_x - A2_rho_x - AB_rho_y - (eps_t./(dt^2) - 0.5./dt).*AA.fv(:,:,1) + 2.*eps_t.*AA.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    % Enforce BC's
    AA.fv(:,1,3) = AA.BC.W;
    AA.fv(end,:,3) = AA.BC.N;    
    
    % Calculate new B
%     [AB_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 2, dr, AB_R.BC, 1);
%     [B2_rho_y, ~, ~] = grad_f((BB.fv(:,:,2).^2)./RHO.fv(:,:,2), 1, dT, BB.BC, 1);
    B2_rho_y = 2.*BB.fv(:,:,2).*BT./RHO.fv(:,:,2) - (BB.fv(:,:,2).^2 .* RHOT)./(RHO.fv(:,:,2).^2);
    AB_rho_x = (AA.fv(:,:,2).*Br./RHO.fv(:,:,2)) + (BB.fv(:,:,2).*Ar)./(RHO.fv(:,:,2)) - (AA.fv(:,:,2).*BB.fv(:,:,2).*RHOr)./(RHO.fv(:,:,2).^2);
    laplace_B = laplace_f(BB.fv(:,:,2), dr, dT, BB.BC, 1);
    [P_y, ~, ~] = grad_f(PP.fv, 1, dT, PP.BC, 1);
    
    BB.fv(:,:,3) = (eps_s.*laplace_B - P_y - AB_rho_x - B2_rho_y - (eps_t./(dt^2) - 0.5./dt).*BB.fv(:,:,1) + 2.*eps_t.*BB.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    BB.fv(:,1,3) = BB.BC.W;
    BB.fv(end,:,3) = BB.BC.N;
    
    % Calculate Residuals
    R_err = (abs(RHO.fv(:,:,3) - RHO.fv(:,:,2)));
    A_err = (abs(AA.fv(:,:,3) - AA.fv(:,:,2)));
    B_err = (abs(BB.fv(:,:,3) - BB.fv(:,:,2)));
    
    if (size(res,1) == 1) && all(res(end,:) == 1)
        res(1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))]; 
    else
        res(end+1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))];
    end
    
    if (size(res, 1) > 500) && (mod(size(res, 1), 2000) == 0)
        fprintf('Iteration Ct: %i\n', size(res, 1));
        fprintf('Current Residual: %0.5e\n', max(res(end, :)));
        toc;
        figure(1);semilogy(1:size(res,1), res(:,1));
        hold on;
        semilogy(1:size(res,1), res(:,2));
        semilogy(1:size(res,1), res(:,3));
        hold off;
        legend('Density', '\rho u', '\rho v');
        fprintf('\n');
    end
    
    if ~isreal(RHO.fv(:,:,end)) || any(any(isnan(RHO.fv(:,:,end))))
        fprintf('Density exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    if ~isreal(AA.fv(:,:,end)) || any(any(isnan(AA.fv(:,:,end))))
        fprintf('Vx exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    if ~isreal(BB.fv(:,:,end)) || any(any(isnan(BB.fv(:,:,end))))
        fprintf('Vy exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
end

%% Post Process

Ux = (AA.fv(:,:,2)./RHO.fv(:,:,2));
Vy = (BB.fv(:,:,2)./RHO.fv(:,:,2));

q2_ij = (Ux).^2 + (Vy).^2;

folderName = ['M_' num2str(M0)];
geomName = [case_name num2str(100*rem(tau,1))];

if ~exist([pwd '\' geomName '\' folderName], 'dir')
    mkdir([pwd '\' geomName '\' folderName]);
end

close all;
figure(1);semilogy(1:size(res,1), res(:,1));
hold on;
semilogy(1:size(res,1), res(:,2));
semilogy(1:size(res,1), res(:,3));
hold off;
legend('Density', '\rho u', '\rho v');
title(['Residual Plot, M=' num2str(M0)]);
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot']);

% plot density
figure();contourf(XX,YY,round(RHO.fv(:,:,2),3), 50)
title(['Density (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\density.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\density']);

figure(); % cp plots
contourf(XX, YY, 1-q2_ij, 50); %./((YY.*cos(XX)).^2)
title(['Pressure Coefficient Contours, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour']);

figure(); % pressure
contourf(XX, YY, PP.fv, 50);
title(['Pressure (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\pressure.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\pressure']);

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
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf']);

% figure(); % field potential
% contourf(XX, YY, round(PHI(:,:,end),5), 50);
% title(['Field Potential, \Phi');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot.pdf']);
% saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(XX, YY, Ux, 50); %./((YY.*cos(XX)).^2)
title(['U velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(XX, YY, Vy, 50); %./((YY.*cos(XX)).^2)
title(['V velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_radius.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_radius']);

% figure();
% contourf(XX, YY, sqrt(M2_ij), 50);
% title(['Mach Number');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\' geomName '\' folderName '\mach.pdf']);
% saveas(gcf, [pwd '\' geomName '\' folderName '\mach']);

% Save Results
save([pwd '\' geomName '\' folderName '\results.mat']);
