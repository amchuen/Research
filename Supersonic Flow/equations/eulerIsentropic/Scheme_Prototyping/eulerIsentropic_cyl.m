clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.005; % spatial diffusion term
eps_t = 0.0013; % time diffusion term
tol = 1e-3;
dt = 0.01^2;
iter_min = 300;
CFL_on = 1;

case_name = 'cylinder';

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.35;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

% Define Grid
dT = 0.01*pi;
dr = 0.01;

r_cyl = 0.5;

% Field Axis Values - body fitted grid
R_range=[   r_cyl+0.5*dr,... % r_cyl
            15];
T_range=[   0,...
            pi];

%% INITIALIZE

GR = struct('R_vals',[], 'T_vals',[], 'RR',[], 'TT',[]);
RHO = struct('fv',[],'BC',[], 'GR', []);
AA = RHO; % radial velocity
BB = RHO; % angular velocity
AB_R = RHO;
PP = RHO;
for i = 1:size(T_range,1)
    GR(i).r_vals = R_range(i,1):dr:R_range(i,2);
    GR(i).T_vals = T_range(i,1):dT:T_range(i,2);
    [GR(i).TT, GR(i).RR] = meshgrid(GR(i).T_vals, GR(i).r_vals);
    GR(i).XX = GR(i).RR .* cos(GR(i).TT);
    GR(i).YY = GR(i).RR .* sin(GR(i).TT);
        
    RHO(i).fv = ones([size(GR(i).TT), 3]); % density
    AA(i).fv = repmat(cos(GR(i).TT).*(1 - (r_cyl^2)./(GR(i).RR.^2)), 1, 1, 3);
    BB(i).fv = repmat(-(1 + (r_cyl^2)./(GR(i).RR.^2)).*sin(GR(i).TT), 1, 1, 3); % rho * v
%     AA(i).fv = RHO(i).fv;
%     BB(i).fv = RHO(i).fv;
    AB_R(i).fv = AA(i).fv .* BB(i).fv ./ RHO(i).fv;
    PP(i).fv = RHO(i).fv(:,:,2).^(gam) ./ (gam .* M0.^2); % pressure
end

%% Boundary Conditions

BC_rho = {'N', 'N', 'D', 'N'}; % west symmetry, east symmetry, far-field, symmetry/wall
% BC_BB = {'N', 'N', 'D', 'N'}; % angular velocity
BC_BB = {'D', 'D', 'D', 'N'}; % angular velocity
BC_AA = {'N', 'N', 'D', 'D'}; % radial velocity
BC_AB_R = {'N', 'N', 'D', 'D'}; % vorticity defined dynamically?

% Continuity
BCval_R = {0, 0, 1, 0};
[RHO.BC, ~] = init_fv(BC_rho, BCval_R, GR.XX);

BCval_rA = {0, 0, (cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2))).*(GR.RR(end,:)+dr), 0};
[rA.BC, ~] = init_fv(BC_AA, BCval_rA, GR.XX);

BCval_BB = {-BB.fv(:,2,2), -BB.fv(:,end-1,2),  -(1 + (r_cyl^2)./((GR.RR(end,:)+dr).^2)).*sin(GR.TT(end,:)), 0};
[BB.BC, B2.BC] = init_fv(BC_BB, BCval_BB, GR.XX);

% Radial momentum
BCval_rA2_R = {0, 0, ((cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2))).^2).*(GR.RR(end,:)+dr), 0};
[rA2_R.BC, ~] = init_fv(BC_AA, BCval_rA2_R, GR.XX);

BCval_AB_R = {0, 0, BCval_BB{3} .* (cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2))), 0};
[AB_R.BC, ~] = init_fv(BC_AB_R, BCval_AB_R, GR.XX); % vorticity term in Radial Momentum

BCval_VV = {0, 0, BCval_BB{3}, 0};
[VV.BC, ~] = init_fv(BC_BB, BCval_VV, GR.XX);

BCval_AA = {0, 0, cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2)), 0};
[AA.BC, ~] = init_fv(BC_AA, BCval_AA, GR.XX);

% Tangential Momentum
BCval_rAB_R = {0, 0, BCval_AB_R{3} .* (GR.RR(end,:)+dr), 0};
[rAB_R.BC, ~] = init_fv(BC_AB_R, BCval_rAB_R, GR.XX); % vorticity term in Angular Momentum

val_test = [0, 0, 1, 0]; 
[PP.BC, ~] = init_fv(BC_rho, num2cell((val_test.^(gam))./(gam .* M0.^2)), GR.XX);

%% START SOLUTION

res = [1, 1, 1];
tic;
while (max(res(end, :)) > tol*max(res(:)))|| (size(res,1) < iter_min) % iterate through time
    
   % Update Field Values
    RHO.fv(:,:,1:2) = RHO.fv(:,:,2:3);
    AA.fv(:,:,1:2) = AA.fv(:,:,2:3);
    BB.fv(:,:,1:2) = BB.fv(:,:,2:3);
    PP.fv = RHO.fv(:,:,2).^(gam) ./ (gam .* M0.^2); % pressure
    
    % Update Boundary Conditions
    BCval_BB = {-BB.fv(:,2,2), -BB.fv(:,end-1,2),  -(1 + (r_cyl^2)./((GR.RR(end,:)+dr).^2)).*sin(GR.TT(end,:)), 0}; % update along axis of symmetry
    [BB.BC, B2.BC] = init_fv(BC_BB, BCval_BB, GR.XX);
    
    BCval_VV = {BCval_BB{1}./RHO.fv(:,2,2), BCval_BB{2}./RHO.fv(:,end-1,2), BCval_BB{3}, 0}; % update along axis of symmetry
    [VV.BC, ~] = init_fv(BC_BB, BCval_VV, GR.XX);

    % Check CFL conditions
    Ur = (AA.fv(:,:,2)./RHO.fv(:,:,2));
    VT = (BB.fv(:,:,2)./RHO.fv(:,:,2))./GR.RR;
    CFL_i = (max(abs(Ur(:)))./(dr) + max(abs(VT(:)))./dT)*dt;

    if CFL_i >= 1.0 && CFL_on
       fprintf('CFL condition not met!\n');
       fprintf('Decreasing time steps!\n');
       dt = dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
    end
    
    % Calculate new RHO
%     [A_r, ~, ~] = grad_f(GR.RR.*AA.fv(:,:,2), 1, dr, AA.BC, 1); % this might be the error
    [A_r, A_rf, A_rb] = grad_f(GR.RR.*AA.fv(:,:,2), 1, dr, rA.BC, 1); 
    [B_T, B_Tf, B_Tb] = grad_f(BB.fv(:,:,2), 2, dT, BB.BC, 1);
	laplace_RHO = laplace_f(RHO.fv(:,:,2), dT, dr, RHO.BC, GR.RR); % note -> need to figure out how to modify radial term for polar coordinates
    RHO.fv(:,:,3) = (eps_s.*laplace_RHO - A_r./GR.RR - B_T./GR.RR - (eps_t./(dt^2) - 0.5./dt).*RHO.fv(:,:,1) + 2.*eps_t.*RHO.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    if (~isreal(RHO.fv(:,:,end)) || any(any(isnan(RHO.fv(:,:,end))))) || any(any(RHO.fv(:,:,end)<0))
        fprintf('Density exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    % Calculate new A
    [rA2_rho_R, ~, ~] = grad_f((GR.RR .* AA.fv(:,:,2).^2)./RHO.fv(:,:,2), 1, dr, rA2_R.BC, 1);
    [AB_rho_T, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 2, dT, AB_R.BC, 1);
    [V_T, ~, ~] = grad_f(BB.fv(:,:,2)./RHO.fv(:,:,2), 2, dT, VV.BC, 1);
    laplace_A = laplace_f(AA.fv(:,:,2)./RHO.fv(:,:,2), dT, dr, AA.BC, GR.RR)...
                - (AA.fv(:,:,2)./RHO.fv(:,:,2))./(GR.RR.^2)...
                -2*V_T./(GR.RR.^2);
    [P_r, ~, ~] = grad_f(PP.fv, 1, dr, PP.BC, 1);
    AA.fv(:,:,3) = (eps_s.*laplace_A - P_r + (BB.fv(:,:,2).^2)./(RHO.fv(:,:,2).*GR.RR) - rA2_rho_R./GR.RR - AB_rho_T./GR.RR - (eps_t./(dt^2) - 0.5./dt).*AA.fv(:,:,1) + 2.*eps_t.*AA.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    if ~isreal(AA.fv(:,:,end)) || any(any(isnan(AA.fv(:,:,end))))
        fprintf('Vr exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
        
    % Calculate new B
    [rAB_rho_R, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2).*GR.RR)./RHO.fv(:,:,2), 1, dr, rAB_R.BC, 1);
    [B2_rho_T, ~, ~] = grad_f((BB.fv(:,:,2).^2)./RHO.fv(:,:,2), 2, dT, B2.BC, 1);
    [U_T, ~, ~] = grad_f(AA.fv(:,:,2)./RHO.fv(:,:,2), 2, dT, AA.BC, 1);
    laplace_B = laplace_f(BB.fv(:,:,2)./RHO.fv(:,:,2), dT, dr, VV.BC, GR.RR)...
                - (BB.fv(:,:,2)./RHO.fv(:,:,2))./(GR.RR.^2)...
                +2*U_T./(GR.RR.^2);
    [P_T, ~, ~] = grad_f(PP.fv, 2, dT, PP.BC, 1);
    BB.fv(:,:,3) = (eps_s.*laplace_B - (BB.fv(:,:,2).*AA.fv(:,:,2))./(RHO.fv(:,:,2).*GR.RR) - P_T./GR.RR - rAB_rho_R./GR.RR - B2_rho_T./GR.RR - (eps_t./(dt^2) - 0.5./dt).*BB.fv(:,:,1) + 2.*eps_t.*BB.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    if any(any(max(abs(BB.fv(:,:,end)./RHO.fv(:,:,end)))>10))
       fprintf('Check Tangential veocity!\n');
    end
    
    if (~isreal(BB.fv(:,:,end)) || any(any(isnan(BB.fv(:,:,end)))))
        fprintf('Vx exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
        
    % Enforce BC's
    RHO.fv(end,:,3) = RHO.BC.N;
    AA.fv(end,:,3) = AA.BC.N;  
    BB.fv(end,:,3) = BB.BC.N;
    AA.fv(1,:,3) = AA.BC.S;
    
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
    
end

%% Post Process

Ur = (AA.fv(:,:,2)./RHO.fv(:,:,2));
VT = (BB.fv(:,:,2)./RHO.fv(:,:,2));

q2_ij = (Ur).^2 + (VT).^2;

folderName = ['M_' num2str(M0)];
geomName = [case_name num2str(100*rem(r_cyl,1))];

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
contourf(GR.XX, GR.YY, 1-q2_ij, 50); %./((YY.*cos(XX)).^2)
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
if M0>1 && solve_half
    plot([fliplr(GR.XX(:,end)'), fliplr(GR.XX(1,:))], 1-[fliplr(q2_ij(:,end)'), fliplr(q2_ij(1,:))]); 
else
    plot([fliplr(GR.XX(:,end)'), fliplr(GR.XX(1,:)), GR.XX(:,1)'], 1-[fliplr(q2_ij(:,end)'), fliplr(q2_ij(1,:)), q2_ij(:,1)']);
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
contourf(GR.XX, GR.YY, Ur, 50); %./((YY.*cos(XX)).^2)
title(['U velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta.pdf']);
saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, VT, 50); %./((YY.*cos(XX)).^2)
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
