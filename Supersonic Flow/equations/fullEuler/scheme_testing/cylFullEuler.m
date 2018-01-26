clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.01; % spatial diffusion term
eps_t = 0.005; % time diffusion term
tol = 1e-5;
dt = 0.0125;
iter_min = 300;
CFL_on = 1;

case_name = 'cylFullEuler';

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 1.1;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

% Define Grid
% dT = 0.025*pi;
% dr = 0.05;

dT = pi/180;
dr = 0.133;

r_cyl = 0.5;

% Field Axis Values - body fitted grid
ranges = [  0.5*pi,  pi;...  % theta
            r_cyl+0.5*dr,... % radius
            75];

GR.r_vals = ranges(2,1):dr:ranges(2,2);
GR.T_vals = ranges(1,1):dT:ranges(1,2);
[GR.TT, GR.RR] = meshgrid(GR.T_vals, GR.r_vals);
GR.XX = GR.RR .* cos(GR.TT);
GR.YY = GR.RR .* sin(GR.TT);
GR.dr = dr;
GR.dT = dT;

% [EE, FF, GG, GR] = manVars('init', ranges, [dT,dr], gam, M0, 1);

%% INITIALIZE

EE = struct('fv',[],'GR', []);
FF = EE; 
GG = EE; 

% Matrix Dimensions
% 1:Y, 2:X, 3:vec, 4:t

EE.fv = repmat(cat(3,   ones(size(GR.XX)),... % rho
                        (cos(GR.TT).*(1 - (r_cyl^2)./(GR.RR.^2))),... % rho * u
                        (-(1 + (r_cyl^2)./(GR.RR.^2)).*sin(GR.TT)),...% rho * v
                        ones(size(GR.XX))./(gam*(gam-1)*M0^2) + 0.5.*((cos(GR.TT).*(1 - (r_cyl^2)./(GR.RR.^2))).^2 + (-(1 + (r_cyl^2)./(GR.RR.^2)).*sin(GR.TT)).^2)),... % rho * E...1/(gam * M0^2).*ones(size(GR.XX))),... % P ... eps_mach error here?
                        1, 1, 1, 3); 

% EE.fv(:,:,5,:) = repmat((EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1),1,1,1,3);
EE.grad_Rf = zeros(size(EE.fv(:,:,:,2)));
EE.grad_Rb = EE.grad_Rf;
EE.grad_Tf = EE.grad_Rf;
EE.grad_Tb = EE.grad_Rf;

FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
FF.fr = zeros(size(FF.fv));

GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
GG.ft = zeros(size(GG.fv));

%% Boundary Condition

% Vector E:
BC.W = {    0,  'N';... rho
            0,  'N';... rho*u -> radial velocity
            0,  'D';... rho*v -> angular velocity...on line of symmetry
            0,  'N';...  rho*E 0,  'N';... % P
            };
            
BC.N = {    ones(size(GR.TT(end,:))), 'D';...  rho
            (cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2))), 'D';... .*(GR.RR(end,:)+dr)
            (-(1 + (r_cyl^2)./((GR.RR(end,:)+dr).^2)).*sin(GR.TT(end,:))), 'D';
            ones(size(GR.TT(end,:)))./(gam*(gam-1)*M0^2)+0.5,  'D';...  rho*E...ones(size(GR.TT(end,:)))./(gam*M0^2), 'D';... % P
%             ones(size(GR.TT(end,:)))./(gam*(gam-1)*M0^2)+0.5.*((cos(GR.TT(end,:)).*(1 - (r_cyl^2)./((GR.RR(end,:)+dr).^2))).^2 + (-(1 + (r_cyl^2)./((GR.RR(end,:)+dr).^2)).*sin(GR.TT(end,:))).^2),  'D';...  rho*E...ones(size(GR.TT(end,:)))./(gam*M0^2), 'D';... % P
            };
            
BC.E = BC.W;
            
BC.S = {    0,  'N';... 
            0,  'D';... % wall is normal to radial direction
            0,  'N';... 
            0,  'N';... 0,  'N',...
            };

%% Solve

res = [];
tic;
rr_n = [0.5.*(GR.RR(2:end,:) + GR.RR(1:(end-1),:)); 0.5.*((GR.RR(end,:)+dr) + GR.RR(end,:))];
rr_s = 0.5.*([2.*GR.RR(1,:)-dr; GR.RR(2:end,:) + GR.RR(1:(end-1),:)]);
alpha = cat(3, zeros(size(GR.RR)), 1./(GR.RR.^2), 1./(GR.RR.^2), zeros(size(GR.RR)));
DF = ((rr_n + rr_s)./(dr^2 .* GR.RR) + 2./(GR.RR.^2 .* dT^2));... + alpha);
while isempty(res) || ((max(res(end, :)) > tol*max(res(res<1)))|| (size(res,1) < iter_min)) % iterate through time
    
    %% Update Field Values
    EE.fv(:,:,:,1:2) = EE.fv(:,:,:,2:3);
%     EE.fv(:,:,5,2) = (EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1);
    PP.fv = (EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1);
    FF.fv = GR.RR.*(EE.fv(:,:,2,2)./EE.fv(:,:,1,2) .* (EE.fv(:,:,:,2)+ reshape([0,0,0,1],1,1,4).*PP.fv));
    GG.fv = (EE.fv(:,:,3,2)./EE.fv(:,:,1,2) .* (EE.fv(:,:,:,2) + reshape([0,0,0,1],1,1,4).*PP.fv));
    MM.fv = cat(3,  zeros(size(EE.fv(:,:,1,2))),...
                    -EE.fv(:,:,3,2).^2 ./ (EE.fv(:,:,1,2).*GR.RR),...   %Rho V^2
                    (EE.fv(:,:,2,2).*EE.fv(:,:,3,2)) ./ (EE.fv(:,:,1,2).*GR.RR),...
                    zeros(size(GR.RR)));   % Rho U V
                
    % Check CFL conditions
    Ur = (EE.fv(:,:,2,2)./EE.fv(:,:,1,2));
    VT = (EE.fv(:,:,3,2)./EE.fv(:,:,1,2));%./GR.RR;
    CFL_i = (max(abs(Ur(:))) + max(abs(VT(:))))*dt./eps_s;
    alpha2 = 2.*dt.*eps_s./(GR.RR.*dT).^2;
    alpha1 = 2.*dt.*eps_s./(GR.RR.*dr.^2).*0.5.*(rr_n+rr_s);
    alphaTot = alpha1 + alpha2;

    if CFL_i ~= 1.0% || (max(abs(alphaTot(:))) > 1)
%        fprintf('CFL condition not met!\n');
%        fprintf('Decreasing time steps!\n');
        cflFactor = 1./CFL_i; %min(1./CFL_i, 1/max(abs(alphaTot(:))));
        dt = dt.*cflFactor;
%        fprintf('New time step:%0.5f\n', dt);
       
        if abs(log10(cflFactor)) > 0.005
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps!\n');
            fprintf('New time step:%0.5e\n', dt);
        end
       
    end
    
    %% Calculate Radial Derivative
    FF.fr(2:end-1,:,:) = (FF.fv(3:end,:,:) - FF.fv(1:end-2,:,:))./(2*GR.dr); % central difference
    FF.fr(end,:,:) = ((GR.RR(end,:)+dr).*(reshape(cell2mat(BC.N(:,1))',[size(FF.fv(end,:,:))]) + reshape([0,0,0,1./(gam.*M0^2)],1,1,4)).*(BC.N{2,1}) - FF.fv(end-1,:,:))./(2*GR.dr);
%     FF.fr(1,:,:) = (FF.fv(:,2,:) - BC.W{2,1}./BC.W{1,1}.*repmat(reshape(cell2mat(BC.W(:,1)),1,1,5),size(FF.fx,1),1,1))./(2*GR.dr);
    FF.fr(1,:,:) = (FF.fv(2,:,:) + FF.fv(1,:,:))./(2*GR.dr);
    
    %% Calculate T Derivative
    GG.ft(:,2:end-1,:) = (GG.fv(:,3:end,:) - GG.fv(:,1:end-2,:))./(2*GR.dT);
    GG.ft(:,end,:) = 0.5.*((GG.fv(:,end,:) - GG.fv(:,end-1,:))./GR.dT + cat(3, (-GG.fv(:,end-1,1) - GG.fv(:,end,1))./GR.dT,...
                                                                            (-GG.fv(:,end-1,2) - GG.fv(:,end,2))./GR.dT,...
                                                                            (GG.fv(:,end-1,3) - GG.fv(:,end,3))./GR.dT,...
                                                                            (-GG.fv(:,end-1,4) - GG.fv(:,end,4))./GR.dT)); % multiply v (B/rho) to BC before operating on boundary
    GG.ft(:,1,:) = 0.5.*((GG.fv(:,2,:) - GG.fv(:,1,:))./GR.dT + (GG.fv(:,1,:) - (4/3.*GG.fv(:,1,:)-1/3.*GG.fv(:,2,:)))./GR.dT);

    
    %% Calculate Laplacians 
    EE.grad_Rf = cat(1, diff(EE.fv(:,:,:,2),1,1)./dr, (reshape(cell2mat(BC.N(:,1))',[size(FF.fv(end,:,:))]) - EE.fv(end,:,:,2))./dr);
    EE.grad_Rb = cat(1, EE.fv(1,:,:,2).*(1 - cat(3, ones(size(GR.RR(1,:))),...
                    -ones(size(GR.RR(1,:))),...
                    (2./(1 - dr./(2*r_cyl)) - 1).*ones(size(GR.RR(1,:))),...
                    ones(size(GR.RR(1,:)))))./(dr), diff(EE.fv(:,:,:,2),1,1)./dr);
    EE.laplace_R = (rr_n.*EE.grad_Rf - rr_s.*EE.grad_Rb)./(GR.dr.*GR.RR);
    
    EE.grad_Tf = cat(2, diff(EE.fv(:,:,:,2),1,2)./dT, (repmat(reshape([1,1,-1,1],1,1,4), size(EE.fv,1), 1, 1).*EE.fv(:,end-1,:,2) - EE.fv(:,end,:,2))./dT);
%     EE.grad_Tb = cat(2, (EE.fv(:,1,:,2) - repmat(reshape([1,1,-1],1,1,3), size(EE.fv,1), 1, 1).*EE.fv(:,2,:,2))./dT, diff(EE.fv(:,:,:,2),1,2)./dT);
%     EE.grad_Tb = cat(2, -3.*EE.fv(:,1,:,2) + 6.*EE.fv(:,2,:,2) -4.*EE.fv(:,3,:,2) + EE.fv(:,4,:,2), diff(EE.fv(:,:,:,2),1,2))./dT;
    EE.grad_Tb = [EE.fv(:,1,:,2) - (4/3.*EE.fv(:,1,:,2) - 1/3.*EE.fv(:,2,:,2)), diff(EE.fv(:,:,:,2),1,2)]./dT;
    EE.grad_T = 0.5.*(EE.grad_Tf + EE.grad_Tb);
    EE.laplace_T = (EE.grad_Tf - EE.grad_Tb)./(dT .* GR.RR.^2);
    
%     EE.laplace = (EE.fx_f(:,:,1:4) - EE.fx_b(:,:,1:4))./(GR.dx) + (EE.fy_f(:,:,1:4) - EE.fy_b(:,:,1:4))./(GR.dy);
    EE.laplace = EE.laplace_R + EE.laplace_T - cat(3,zeros(size(GR.RR)), EE.fv(:,:,2:3,2), zeros(size(GR.RR)))./(GR.RR.^2) + (2.*cat(3, zeros(size(GR.RR)), -EE.grad_T(:,:,3), EE.grad_T(:,:,2), zeros(size(GR.RR))))./(GR.RR.^2);
    
    %% Pressure Terms
    PP.grad_Rf = [diff(PP.fv, 1, 1)./dr; (1./(gam*M0^2) - PP.fv(end,:))./dr];
%     PP.grad_Rb = [zeros(size(PP.fv(1,:))); diff(PP.fv, 1, 1)./dr];
    PP.grad_Rb = [GG.fv(1,:,3)./r_cyl; diff(PP.fv, 1, 1)./dr];
    PP.grad_R = 0.5*(PP.grad_Rf + PP.grad_Rb);
    
    PP.grad_Tf = [diff(PP.fv, 1,2)./dT, (PP.fv(:,end-1) - PP.fv(:,end))./dT];
%     PP.grad_Tb = [(PP.fv(:,1) - PP.fv(:,2))./dT, diff(PP.fv, 1,2)./dT];
%     PP.grad_Tb = [(3.*PP.fv(:,2) - 2.*PP.fv(:,1) - PP.fv(:,3))./dT, diff(PP.fv, 1,2)./dT];
    PP.grad_Tb = [(PP.fv(:,1) - (4/3.*PP.fv(:,1)-1/3.*PP.fv(:,2)))./dT, diff(PP.fv, 1,2)./dT];
    PP.grad_T = 0.5*(PP.grad_Tf + PP.grad_Tb); 
    
%     PP.grad = cat(3, zeros(size(GR.RR)), PP.grad_R, PP.grad_T./GR.RR);
    
    %% Calculate Next Time Step
    
%     DF = (2./(GR.dx^2) + 2./(GR.dy^2)); % Dufort-Frankel coefficient
    pTerms = cat(3, zeros(size(GR.XX)),...
                        PP.grad_R,...
                        PP.grad_T./GR.RR,...
                        zeros(size(GR.XX)));
%     EE.fv(:,:,:,3) = (eps_s.*EE.laplacian - PP.grad - (FF.grad_R + GG.grad_T)./RR - MM.fv + 0.5.*EE.fv(:,:,:,1)./dt)./...
%                     (0.5./dt + 0.5.*eps_s.*DF_coeff);
    EE.fv(:,:,1:4,3) = (eps_s.*(EE.laplace + DF.*EE.fv(:,:,1:4,2)) - 0.5*EE.fv(:,:,1:4,1).*(eps_s*DF - 1/dt) - (FF.fr(:,:,1:4) + GG.ft(:,:,1:4) + pTerms) - MM.fv)./(0.5/dt + 0.5*eps_s*DF);
    
    %% Calculate Residual
    
    if (size(res,1) == 1) && all(res(end,:) == 0)
        res(1,:) = reshape(max(max(abs(EE.fv(:,:,1:4,3) - EE.fv(:,:,1:4,2)))),1,4,1);
    else
        res(end+1,:) = reshape(max(max(abs(EE.fv(:,:,1:4,3) - EE.fv(:,:,1:4,2)))),1,4,1);
    end

%     figure(1);semilogy(1:size(res,1), res(:,1));
%     hold on;
%     semilogy(1:size(res,1), res(:,2));
%     semilogy(1:size(res,1), res(:,3));
%     semilogy(1:size(res,1), res(:,4));
%     hold off;
%     legend('Density', '\rho u', '\rho v', '\rho\epsilon', 'Location' ,'bestoutside');
    figure(1);contourf(GR.XX, GR.YY, EE.fv(:,:,4,end), 50); colorbar;
    drawnow;
    if (size(res, 1) > 500) && (mod(size(res, 1), 2000) == 0)
        fprintf('Iteration Ct: %i\n', size(res, 1));
        fprintf('Current Residual: %0.5e\n', max(res(end, :)));
        toc;
        fprintf('\n');
    end
    
    if (~isreal(EE.fv(:,:,:,end)) || any(any(any(isnan(EE.fv(:,:,:,end))))))
        error('Solution exhibits non-solutions (either non-real or NaN) in nodes!\n');
%         break;
    end
    
end

%% Post Process

rho = EE.fv(:,:,1,2);
aa = EE.fv(:,:,2,2);
bb = EE.fv(:,:,3,2);
cc = EE.fv(:,:,4,2);
pressure = EE.fv(:,:,end,2);

levels = 50;

Ur = (aa./rho);
VT = (bb./rho);

q2_ij = (Ur).^2 + (VT).^2;

folderName = ['M_' num2str(M0)];
geomName = [case_name num2str(100*rem(tau,1))];
sizeName = [num2str(size(GR.XX,1)) '_' num2str(size(GR.XX,2)) '_' num2str(round(dx,3)) '_' num2str(round(dy,3))];
% dirName = [pwd '\' geomName '\' sizeName '\' folderName];
dirName = [pwd '\' geomName '\' folderName];

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

close all;
for i = 1:size(res,2)
    figure(1);semilogy(1:size(res,1), res(:,i));
    hold on;
end
hold off;
legend('Density', '\rho u', '\rho v', '\rho \epsilon');
title(['Residual Plot, M=' num2str(M0)]);
xlabel('# of iterations');
ylabel('Residual (Error)');

saveas(gcf, [dirName '\residual_plot.pdf']);
saveas(gcf, [dirName '\residual_plot']);

% plot density
figure();contourf(GR.XX,GR.YY,rho, levels)
title(['Density (Normalized), M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [dirName '\density.pdf']);
saveas(gcf, [dirName '\density']);

figure(); % cp plots
contourf(GR.XX, GR.YY, 1-q2_ij, levels); %./((YY.*cos(XX)).^2)
title(['Pressure Coefficient Contours, M=' num2str(M0)]);
colorbar('eastoutside');
axis equal
saveas(gcf, [dirName '\cp_contour.pdf']);
saveas(gcf, [dirName '\cp_contour']);

figure(); % pressure
contourf(GR.XX, GR.YY, pressure, levels);
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
title(['C_p on surface of Cylinder, M=' num2str(M0)]);
set(gca, 'Ydir', 'reverse');
% axis equal
saveas(gcf, [dirName '\cp_surf.pdf']);
saveas(gcf, [dirName '\cp_surf']);

% figure(); % field potential
% contourf(XX, YY, round(PHI(:,:,end),5));
% title(['Field Potential, \Phi');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [dirName '\phi_pot.pdf']);
% saveas(gcf, [dirName '\phi_pot']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, Ur, levels); %./((YY.*cos(XX)).^2)
title(['U velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [dirName '\phi_theta.pdf']);
saveas(gcf, [dirName '\phi_theta']);

figure(); % theta-dir velocity plots
contourf(GR.XX, GR.YY, VT, levels); %./((YY.*cos(XX)).^2)
title(['V velocity, M=' num2str(M0)]);
colorbar('eastoutside');
saveas(gcf, [dirName '\phi_radius.pdf']);
saveas(gcf, [dirName '\phi_radius']);

% figure();
% contourf(XX, YY, sqrt(M2_ij));
% title(['Mach Number');
% colorbar('eastoutside');
% axis equal
% saveas(gcf, [dirName '\mach.pdf']);
% saveas(gcf, [dirName '\mach']);

% Save Results
save([dirName '\results.mat']);

