clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.005; % spatial diffusion term
eps_t = 0.0013; % time diffusion term
tol = 1e-4;
dt = 0.002; %0.01;
iter_min = 300;
CFL_on = 1;

case_name = 'biconvex_test';

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.98;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dx = 0.01;
dy = 0.08/5;

% Field Axis Values - body fitted grid
x_range=[   -7-39*dx,...
            7 + 20*dx];
y_range=[   0,...
            10];
        
%% INITIALIZE

GR = struct('x_vals',[], 'y_vals',[], 'XX',[], 'YY',[]);
RHO = struct('fv',[],'BC',[]);
AA = RHO;
BB = RHO;
AB_R = RHO;
PP = RHO;
for i = 1:size(x_range,1)
    GR(i).x_vals = x_range(i,1):dx:x_range(i,2);
    GR(i).y_vals = y_range(i,1):dy:y_range(i,2);
    [GR(i).XX, GR(i).YY] = meshgrid(GR(i).x_vals, GR(i).y_vals);
        
    RHO(i).fv = ones([size(GR(i).XX), 3]); % density
    AA(i) = RHO(i); % rho * u
    BB(i).fv = zeros(size(RHO(i).fv)); % rho * v
    AB_R(i).fv = AA(i).fv .* BB(i).fv ./ RHO(i).fv;
    PP(i).fv = RHO(i).fv(:,:,2).^(gam) ./ (gam .* M0.^2); % pressure
end

%% Body Values
tau = 0.1;
YY_B = [zeros(size(GR.x_vals(GR.x_vals <0))), ...
        2*tau.*GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1)).*(1- GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1))),...
        zeros(size(GR.x_vals(GR.x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i+1) = (YY_B(i) - YY_B(i-1))/(2*dx);
end

%% Boundary Conditions

BC_rho = {'D', 'N', 'D', 'N'}; % inlet, outlet, far-field, symmetry/wall
BC_AA = {'D', 'N', 'D', 'N'};        
BC_BB = {'D', 'N', 'D', 'D'};
BC_AB_R = {'D', 'N', 'D', 'D'};

val_test = [1, 0, 1, 0]; 
BCval = {1, 0, 1, 0}; 
BCval_BB = {0, 0, 0, RHO.fv(1,:,2).*dyBdx};
BCval_AB_R = {0, 0, 0, AA.fv(1,:, 2).*dyBdx};

RHO.BC = init_fv(BC_rho, BCval, GR.XX);
AA.BC = init_fv(BC_AA, BCval, GR.XX);
BB.BC = init_fv(BC_BB, BCval_BB, GR.XX);
AB_R.BC = init_fv(BC_AB_R, BCval_AB_R, GR.XX);
PP.BC = init_fv(BC_rho, num2cell((val_test.^(gam))./(gam .* M0.^2)), GR.XX);

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
%     BCval_BB = {0, 0, 0, RHO.fv(1,:,2).*dyBdx};
%     BCval_AB_R = {0, 0, 0, AA.fv(1,:, 2).*dyBdx};
%     BB.BC = init_fv(BC_BB, BCval_BB, GR.XX);
%     AB_R.BC = init_fv(BC_AB_R, BCval_AB_R, GR.XX);
    BB.BC.S = RHO.fv(1,:,2).*dyBdx;
    AB_R.BC.S = AA.fv(1,:, 2).*dyBdx;
    
    % Check CFL conditions
    Ux = (AA.fv(:,:,2)./RHO.fv(:,:,2));
    Vy = (BB.fv(:,:,2)./RHO.fv(:,:,2));
    CFL_i = (max(abs(Ux(:)))./dx + max(abs(Vy(:)))./dy)*dt;

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
    
    % Calculate new RHO
    [Ax, ~, ~] = grad_f(AA.fv(:,:,2), 2, dx, AA.BC, 1);
    [By, ~, ~] = grad_f(BB.fv(:,:,2), 1, dy, BB.BC, 1);
	laplace_RHO = laplace_f(RHO.fv(:,:,2), dx, dy, RHO.BC, 1);
    RHO.fv(:,:,3) = (eps_s.*laplace_RHO - Ax - By - (eps_t./(dt^2) - 0.5./dt).*RHO.fv(:,:,1) + 2.*eps_t.*RHO.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
        
    % Calculate new A
    [A2_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).^2)./RHO.fv(:,:,2), 2, dx, AA.BC, 1);
    [AB_rho_y, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 1, dy, AB_R.BC, 1);
    laplace_A = laplace_f(AA.fv(:,:,2), dx, dy, AA.BC, 1);
    [P_x, ~, ~] = grad_f(PP.fv, 2, dx, PP.BC, 1);
    AA.fv(:,:,3) = (eps_s.*laplace_A - P_x - A2_rho_x - AB_rho_y - (eps_t./(dt^2) - 0.5./dt).*AA.fv(:,:,1) + 2.*eps_t.*AA.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    % Calculate new B
    [AB_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 2, dx, AB_R.BC, 1);
    [B2_rho_y, ~, ~] = grad_f((BB.fv(:,:,2).^2)./RHO.fv(:,:,2), 1, dy, BB.BC, 1);
    laplace_B = laplace_f(BB.fv(:,:,2), dx, dy, BB.BC, 1);
    [P_y, ~, ~] = grad_f(PP.fv, 1, dy, PP.BC, 1);
    BB.fv(:,:,3) = (eps_s.*laplace_B - P_y - AB_rho_x - B2_rho_y - (eps_t./(dt^2) - 0.5./dt).*BB.fv(:,:,1) + 2.*eps_t.*BB.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    % Enforce Dirichlet BC's
    RHO.fv(:,1,3) = RHO.BC.W;
    RHO.fv(end,:,3) = RHO.BC.N;
    AA.fv(:,1,3) = AA.BC.W;
    AA.fv(end,:,3) = AA.BC.N; 
    BB.fv(:,1,3) = BB.BC.W;
    BB.fv(end,:,3) = BB.BC.N;
    BB.fv(1,:,3) = BB.BC.S;
    
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
figure();contourf(GR.XX,GR.YY,round(RHO.fv(:,:,2),3), 50)
title(['Density (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\density.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\density']);

figure(); % cp plots
contourf(GR.XX, GR.YY, 1-q2_ij, 50); %./((RR.*cos(TT)).^2)
title(['Pressure Coefficient Contours, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour']);

figure(); % pressure
contourf(GR.XX, GR.YY, PP.fv, 50);
title(['Pressure (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [pwd '\' geomName '\' folderName '\pressure.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\pressure']);

figure();
% plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
plot(GR.XX(1,:), 1-q2_ij(1,:));
xlabel('X');
%     ylabel('\phi_{\theta}');
ylabel('C_p');
title(['C_p on surface, M=' num2str(M0)]);
set(gca, 'Ydir', 'reverse');
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf']);

% figure(); % field potential
% contourf(GR.XX, GR.YY, round(PHI(:,:,end),5), 50);
% title(['Field Potential, \Phi');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot.pdf']);
% saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, Ux, 50); %./((RR.*cos(TT)).^2)
title(['U velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, Vy, 50); %./((RR.*cos(TT)).^2)
title(['V velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_radius.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_radius']);

% figure();
% contourf(GR.XX, GR.YY, sqrt(M2_ij), 50);
% title(['Mach Number');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [pwd '\' geomName '\' folderName '\mach.pdf']);
% saveas(gcf, [pwd '\' geomName '\' folderName '\mach']);

% Save Results
save([pwd '\' geomName '\' folderName '\results.mat']);

%% Save copy of file

copyfile('biconvex.m',[pwd '\' geomName '\' folderName '\biconvex.m']);