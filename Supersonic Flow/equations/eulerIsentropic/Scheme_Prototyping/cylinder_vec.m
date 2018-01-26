clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.0725; % spatial diffusion term
eps_t = 0.005; % time diffusion term
tol = 1e-4;
dt = 0.1;
iter_min = 300;
CFL_on = 1;

case_name = 'cylinder_vectorized';

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 0.53;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

% Define Grid
dT = 0.025*pi;
dr = 0.05;

r_cyl = 0.5;

% Field Axis Values - body fitted grid
R_range=[   r_cyl+0.5*dr,... % r_cyl
            15];
T_range=[   0,...
            pi];
        
%% INITIALIZE

GR = struct('R_vals',[], 'T_vals',[], 'RR',[], 'TT',[]);
EE = struct('fv',[],'BC',[], 'GR', []);
FF = EE; % radial velocity
GG = EE; % angular velocity

GR.r_vals = R_range(1):dr:R_range(2);
GR.T_vals = T_range(1):dT:T_range(2);
[GR.TT, GR.RR] = meshgrid(GR.T_vals, GR.r_vals);
GR.XX = GR.RR .* cos(GR.TT);
GR.YY = GR.RR .* sin(GR.TT);

EE.fv = repmat(cat(3,  ones(size(GR.TT)),... % density
            (cos(GR.TT).*(1 - (r_cyl^2)./(GR.RR.^2))),... % rho * u
            (-(1 + (r_cyl^2)./(GR.RR.^2)).*sin(GR.TT))), 1, 1, 1, 3); % rho * v
FF.fv = cat(3,  GR.RR.*EE.fv(:,:,2,2),...   % rA
                GR.RR.*(EE.fv(:,:,2,2).^2)./EE.fv(:,:,1,2),...  % rA^2/Rho
                GR.RR.*EE.fv(:,:,2,2).*EE.fv(:,:,3,2)./EE.fv(:,:,1,2)); % rAB/Rho
GG.fv = cat(3,  EE.fv(:,:,3,2),...  % B
                (EE.fv(:,:,2,2).*EE.fv(:,:,3,2))./EE.fv(:,:,1,2),...    % AB/Rho
                (EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2));   % B^2/Rho
MM.fv = cat(3,  zeros(size(EE.fv(:,:,1,2))),...
                -EE.fv(:,:,3,2).^2 ./ (EE.fv(:,:,1,2).*GR.RR),...   %Rho V^2
                (EE.fv(:,:,2,2).*EE.fv(:,:,3,2)) ./ (EE.fv(:,:,1,2).*GR.RR));   % Rho U V
PP.fv = EE.fv(:,:,1,2).^(gam) ./ (gam .* M0.^2); % pressure

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

BCval_BB = {0, 0,  -(1 + (r_cyl^2)./((GR.RR(end,:)+dr).^2)).*sin(GR.TT(end,:)), 0};
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

res = [0, 0, 0];
tic;
rr_n = repmat([0.5.*(GR.RR(2:end,:) + GR.RR(1:(end-1),:)); 0.5.*((GR.RR(end,:)+dr) + GR.RR(end,:))],1,1,3);
rr_s = repmat(0.5.*([2.*GR.RR(1,:)-dr; GR.RR(2:end,:) + GR.RR(1:(end-1),:)]),1,1,3);
RR = repmat(GR.RR, 1,1,3);
while (max(res(end, :)) > tol*max(res(res<1)))|| (size(res,1) < iter_min) % iterate through time
    
    %% Update
    EE.fv(:,:,:,1:2) = EE.fv(:,:,:,2:3);
    FF.fv = cat(3,  GR.RR.*EE.fv(:,:,2,2),...   % rA
                    GR.RR.*(EE.fv(:,:,2,2).^2)./EE.fv(:,:,1,2),...  % rA^2/Rho
                    GR.RR.*EE.fv(:,:,2,2).*EE.fv(:,:,3,2)./EE.fv(:,:,1,2)); % rAB/Rho
    GG.fv = cat(3,  EE.fv(:,:,3,2),...  % B
                    (EE.fv(:,:,2,2).*EE.fv(:,:,3,2))./EE.fv(:,:,1,2),...    % AB/Rho
                    (EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2));   % B^2/Rho
    MM.fv = cat(3,  zeros(size(EE.fv(:,:,1,2))),...
                    -EE.fv(:,:,3,2).^2 ./ (EE.fv(:,:,1,2).*GR.RR),...   %Rho V^2
                    (EE.fv(:,:,2,2).*EE.fv(:,:,3,2)) ./ (EE.fv(:,:,1,2).*GR.RR));   % Rho U V
    PP.fv = EE.fv(:,:,1,2).^(gam) ./ (gam .* M0.^2); % pressure
    
    % Check CFL conditions
    Ur = (EE.fv(:,:,2,2)./EE.fv(:,:,1,2));
    VT = (EE.fv(:,:,3,2)./EE.fv(:,:,1,2))./GR.RR;
    CFL_i = (max(abs(Ur(:)))./(dr) + max(abs(VT(:)))./dT)*dt;

    if CFL_i >= 1.0 && CFL_on
       fprintf('CFL condition not met!\n');
       fprintf('Decreasing time steps!\n');
       dt = dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
    end
    
    %% Calculate Radial Derivative
    FF.grad_Rf = cat(1, diff(FF.fv,1,1)./dr,(cat(3, rA.BC.N, rA2_R.BC.N, rAB_R.BC.N)  - FF.fv(end,:,:))./dr);
    FF.grad_Rb = cat(1, (FF.fv(1,:,:) - cat(3, rA.BC.S, rA2_R.BC.S, rAB_R.BC.S))./(0.5*dr), diff(FF.fv,1,1)./dr);
    FF.grad_R = 0.5.*(FF.grad_Rf + FF.grad_Rb);
    
    %% Calculate Angular Derivatives
    GG.BC.Tf = cat(3, (-GG.fv(:,end-1,1) - GG.fv(:,end,1))./dT,...
                (-GG.fv(:,end-1,2) - GG.fv(:,end,2))./dT,...
                (GG.fv(:,end-1,3) - GG.fv(:,end,3))./dT);
    GG.BC.Tb = cat(3, (GG.fv(:,1,1) + GG.fv(:,2,1))./dT,...
                (GG.fv(:,1,2) + GG.fv(:,2,2))./dT,...
                (GG.fv(:,1,3) - GG.fv(:,2,3))./dT);
    GG.grad_T = 0.5.*(cat(2, diff(GG.fv,1,2)./dT, GG.BC.Tf) + cat(2, GG.BC.Tb, diff(GG.fv,1,2)./dT));
    
    %% Calculate Laplacian
    EE.grad_Rf = cat(1, diff(EE.fv(:,:,:,2),1,1)./dr, (cat(3, RHO.BC.N, AA.BC.N, BB.BC.N) - EE.fv(end,:,:,2))./dr);
    EE.grad_Rb = cat(1, EE.fv(1,:,:,2).*(1 - cat(3, ones(size(GR.RR(1,:))),...
                    (0.5*dr - 2*r_cyl)./(0.5*dr + 2*r_cyl).*((GR.RR(1,:))./(GR.RR(1,:)-dr)).^2,...
                    (GR.RR(1,:))./(GR.RR(1,:)-dr)))./(dr), diff(EE.fv(:,:,:,2),1,1)./dr);
    EE.laplace_R = (rr_n.*EE.grad_Rf - rr_s.*EE.grad_Rb)./(dr.*RR);
    
    EE.grad_Tf = cat(2, diff(EE.fv(:,:,:,2),1,2)./dT, (repmat(reshape([1,1,-1],1,1,3), size(EE.fv,1), 1, 1).*EE.fv(:,end-1,:,2) - EE.fv(:,end,:,2))./dT);
    EE.grad_Tb = cat(2, (EE.fv(:,1,:,2) - repmat(reshape([1,1,-1],1,1,3), size(EE.fv,1), 1, 1).*EE.fv(:,2,:,2))./dT, diff(EE.fv(:,:,:,2),1,2)./dT);
    EE.grad_T = 0.5.*(EE.grad_Tf + EE.grad_Tb);
    EE.laplace_T = (EE.grad_Tf - EE.grad_Tb)./(dT .* RR.^2);
    
%     EE.laplacian = EE.laplace_R + EE.laplace_T + ((rr_n + rr_s)./(dr^2 .* RR) + 2./(RR.^2 .* dT^2)).*(EE.fv(:,:,:,2)-0.5.*EE.fv(:,:,:,1))... % scalar laplacian
%                    - cat(3,zeros(size(GR.RR)), 0.5.*EE.fv(:,:,2:3,1))./(RR.^2)... % E-tilde
%                    + (2.*cat(3, zeros(size(GR.RR)), -EE.grad_T(:,:,3), EE.grad_T(:,:,2)))./(RR.^2); % angular deriv for vector laplacian
    EE.laplacian = EE.laplace_R + EE.laplace_T - cat(3,zeros(size(GR.RR)), EE.fv(:,:,2:3,2))./(RR.^2) + (2.*cat(3, zeros(size(GR.RR)), -EE.grad_T(:,:,3), EE.grad_T(:,:,2)))./(RR.^2);
    
    %% Pressure Terms
    PP.grad_Rf = [diff(PP.fv, 1, 1)./dr; (PP.BC.N - PP.fv(end,:))./dr];
    PP.grad_Rb = [zeros(size(PP.fv(1,:))); diff(PP.fv, 1, 1)./dr];
    PP.grad_R = 0.5*(PP.grad_Rf + PP.grad_Rb);
    
    PP.grad_Tf = [diff(PP.fv, 1,2)./dT, (PP.fv(:,end-1) - PP.fv(:,end))./dT];
    PP.grad_Tb = [(PP.fv(:,1) - PP.fv(:,2))./dT, diff(PP.fv, 1,2)./dT];
    PP.grad_T = 0.5*(PP.grad_Tf + PP.grad_Tb); 
    
    PP.grad = cat(3, zeros(size(GR.RR)), PP.grad_R, PP.grad_T./GR.RR);
    
    %% Step forward in time
%     EE.fv(:,:,:,3) = (eps_s .* EE.laplacian - PP.grad - (FF.grad_R + GG.grad_T)./RR - MM.fv)./...
%                     (0.5/dt + 0.5.*eps_s.*((rr_n + rr_s)./(RR.*dr^2)+2./(RR.^2 .* dT^2) + cat(3,zeros(size(GR.RR)), 1./RR(:,:,2:3))));
    
    EE.fv(:,:,:,3) = (eps_s .* EE.laplacian - PP.grad - (FF.grad_R + GG.grad_T)./RR - MM.fv + 2*eps_t./(dt^2).*EE.fv(:,:,:,2) - (eps_t./(dt^2) - 0.5/dt).*EE.fv(:,:,:,1))./(eps_t./(dt^2) + 0.5/dt);
    
    if any(any(max(abs(EE.fv(:,:,3,end)./EE.fv(:,:,1,end)))>10))
       fprintf('Check Tangential veocity!\n');
    end
    
    if (~isreal(EE.fv(:,:,:,end)) || any(any(any(isnan(EE.fv(:,:,:,end))))))
        fprintf('Solution exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    %% Enforce BC's
%     EE.fv(end,:,1,3) = RHO.BC.N;
%     EE.fv(end,:,2,3) = AA.BC.N;  
%     EE.fv(end,:,3,3) = BB.BC.N;
%     EE.fv(:,1,3,3) = zeros(size(EE.fv(:,1,3,3)));
%     EE.fv(:,end,3,3) = zeros(size(EE.fv(:,end,3,3)));
%     AA.fv(1,:,3) = AA.BC.S;
    
    % Calculate Residuals
    R_err = (abs(EE.fv(:,:,1,3) - EE.fv(:,:,1,2)));
    A_err = (abs(EE.fv(:,:,2,3) - EE.fv(:,:,2,2)));
    B_err = (abs(EE.fv(:,:,3,3) - EE.fv(:,:,3,2)));
    
    if (size(res,1) == 1) && all(res(end,:) == 0)
        res(1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))]; 
    else
        res(end+1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))];
    end
    figure(1);contourf(GR.XX, GR.YY, EE.fv(:,:,1,end), 50);colorbar;axis equal;drawnow;
    if (size(res, 1) > 500) && (mod(size(res, 1), 2000) == 0)
        fprintf('Iteration Ct: %i\n', size(res, 1));
        fprintf('Current Residual: %0.5e\n', max(res(end, :)));
        toc;
%         figure(1);semilogy(1:size(res,1), res(:,1));
%         hold on;
%         semilogy(1:size(res,1), res(:,2));
%         semilogy(1:size(res,1), res(:,3));
%         hold off;
%         legend('Density', '\rho u', '\rho v', 'Location' ,'bestoutside');
%         figure();contourf(GR.XX, GR.YY, EE.fv(:,:,1,end), 50);
        fprintf('\n');
    end
    
end

%% Post Process

rho = EE.fv(:,:,1,2);
aa = EE.fv(:,:,2,2);
bb = EE.fv(:,:,3,2);

Ur = (aa./rho);
VT = (bb./rho);

q2_ij = (Ur).^2 + (VT).^2;

folderName = ['M_' num2str(M0)];
geomName = [case_name num2str(100*rem(r_cyl,1))];
sizeName = [num2str(size(GR.RR,1)) '_' num2str(size(GR.RR,2)) '_' num2str(round(dr,3)) '_' num2str(round(dT,3))];
dirName = [pwd '\' geomName '\' sizeName '\' folderName];

if ~exist(dirName, 'dir')
    mkdir(dirName);
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

saveas(gcf, [dirName '\residual_plot.pdf']);
saveas(gcf, [dirName '\residual_plot']);

% plot density
figure();contourf(GR.XX,GR.YY,rho, 50)
title(['Density (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [dirName '\density.pdf']);
saveas(gcf, [dirName '\density']);

figure(); % cp plots
contourf(GR.XX, GR.YY, 1-q2_ij, 50); %./((YY.*cos(XX)).^2)
title(['Pressure Coefficient Contours, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [dirName '\cp_contour.pdf']);
saveas(gcf, [dirName '\cp_contour']);

figure(); % pressure
contourf(GR.XX, GR.YY, PP.fv, 50);
title(['Pressure (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [dirName '\pressure.pdf']);
saveas(gcf, [dirName '\pressure']);

figure();
solve_half = 0;
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
saveas(gcf, [dirName '\cp_surf.pdf']);
saveas(gcf, [dirName '\cp_surf']);

% figure(); % field potential
% contourf(XX, YY, round(PHI(:,:,end),5), 50);
% title(['Field Potential, \Phi');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [dirName '\phi_pot.pdf']);
% saveas(gcf, [dirName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, Ur, 50); %./((YY.*cos(XX)).^2)
title(['U velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [dirName '\phi_theta.pdf']);
saveas(gcf, [dirName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, VT, 50); %./((YY.*cos(XX)).^2)
title(['V velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [dirName '\phi_radius.pdf']);
saveas(gcf, [dirName '\phi_radius']);

% figure();
% contourf(XX, YY, sqrt(M2_ij), 50);
% title(['Mach Number');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [dirName '\mach.pdf']);
% saveas(gcf, [dirName '\mach']);

% Save Results
save([dirName '\results.mat']);
