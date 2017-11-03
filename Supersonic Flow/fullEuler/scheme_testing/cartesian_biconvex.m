clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 0.005; % spatial diffusion term
eps_t = 0.005; % time diffusion term
tol = 1e-5;
dt = 0.0125;
iter_min = 300;
CFL_on = 1;

case_name = 'biconvexFullEuler';

%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 1.4;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

% Define Grid
dx = 0.05;
dy = 0.08;

% Field Axis Values - body fitted grid
% Field Axis Values - body fitted grid
% x_range=[   -4-39*dx,...
%             7 + 20*dx];
x_range=[   -1,...
            4];
y_range=[   0,...
            10];
        
%% INITIALIZE

GR = struct('Y_vals',[], 'X_vals',[], 'XX',[], 'YY',[]);
EE = struct('fv',[],'GR', []);
FF = EE; 
GG = EE; 

GR.y_vals = y_range(1):dy:y_range(2);
GR.x_vals = x_range(1):dx:x_range(2);
GR.dy = dy; GR.dx = dx;

[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

% Matrix Dimensions
% 1:Y, 2:X, 3:vec, 4:t

EE.fv = repmat(cat(3,   ones(size(GR.XX)),... % rho
                        ones(size(GR.XX)),... % rho * u
                        zeros(size(GR.XX)),...% rho * v
                        ones(size(GR.XX))./(gam*(gam-1)*M0^2) + 0.5,... % rho * E
                        1/(gam * M0^2).*ones(size(GR.XX))),... % P ... eps_mach error here?
                        1, 1, 1, 3); 

% EE.fv(:,:,5,:) = repmat((EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1),1,1,1,3);
EE.fx_f = zeros(size(EE.fv(:,:,:,2)));
EE.fx_b = EE.fx_f;
EE.fy_f = EE.fx_f;
EE.fy_b = EE.fx_f;

FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
FF.fx = zeros(size(FF.fv));

GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
GG.fy = zeros(size(GG.fv));
% MM.fv = zeros(size(FF.fv)); % stores pressure gradients for momentum conservation

%% Body Values
tau = 0.1;
YY_B = [zeros(size(GR.x_vals(GR.x_vals <0))), ...
        2*tau.*GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1)).*(1- GR.x_vals((GR.x_vals>=0)&(GR.x_vals <=1))),...
        zeros(size(GR.x_vals(GR.x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i+1) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

%% Boundary Condition

% Vector E:
BC.W = {    1,  'D';... rho
            1,  'D';... rho*u
            0,  'D';... rho*v
            1./(gam*(gam-1)*M0^2)+0.5,  'D';...  rho*E
            1/(gam*M0^2), 'D';... % P
            };
            
BC.N = BC.W;
if M0 > 1          
    BC.E = {    0,  'N';... 
                0,  'N';...
                0,  'N';...
                0,  'N';...
                0,  'N';...
                };
else
    BC.E = BC.W;
end
            
BC.S = {    0,  'N';... 
            0,  'N';...
            dyBdx, 'D';... % wall is normal to this direction, i.e. wall is parallel to x-dir
            0,  'N';...
            0,  'N',...
            };

%% Solve

res = [];
tic;

while isempty(res) || ((max(res(end, :)) > tol*max(res(res<1)))|| (size(res,1) < iter_min)) % iterate through time
    
    %% Update Field Values
    EE.fv(:,:,:,1:2) = EE.fv(:,:,:,2:3);
    EE.fv(:,:,5,2) = (EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1);
    FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
    GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
                
    % Check CFL conditions
    Ux = (EE.fv(:,:,2,2)./EE.fv(:,:,1,2));
    Vy = (EE.fv(:,:,3,2)./EE.fv(:,:,1,2));
    CFL_i = (max(abs(Ux(:)))./(dx) + max(abs(Vy(:)))./dy)*dt;

    if CFL_i >= 1.0
       fprintf('CFL condition not met!\n');
       if CFL_on
           fprintf('Decreasing time steps!\n');
           dt = dt*0.8 / CFL_i;
           fprintf('New time step:%0.5e\n', dt);
       end
       fprintf('\n');
    end
    
    %% Calculate X Derivative
    FF.fx(:,2:end-1,:) = (FF.fv(:,3:end,:) - FF.fv(:,1:end-2,:))./(2*GR.dx); % central difference
    if M0 > 1
        FF.fx(:,end,:) = repmat(reshape(cell2mat(BC.E(:,1)),1,1,5),size(FF.fx,1),1,1);
    else
        FF.fx(:,end,:) = (BC.E{2,1}./BC.E{1,1}.*repmat(reshape(cell2mat(BC.E(:,1)),1,1,5),size(FF.fx,1),1,1) - FF.fv(:,end-1,:))./(2*GR.dx);
    end
    FF.fx(:,1,:) = (FF.fv(:,2,:) - BC.W{2,1}./BC.W{1,1}.*repmat(reshape(cell2mat(BC.W(:,1)),1,1,5),size(FF.fx,1),1,1))./(2*GR.dx);

    EE.fx_f(:,1:end-1,:) = diff(EE.fv(:,:,:,2),1,2)./GR.dx; % forward difference
    EE.fx_b(:,2:end,:) = EE.fx_f(:,1:end-1,:); % backward difference
    if M0 > 1 % Neumann condition for supersonic flow
        EE.fx_f(:,end,:) = (EE.fv(:,end-1,:,2) - EE.fv(:,end,:,2))./GR.dx; % right boundary
    else % Dirichlet/Free-Stream Condition for Subsonic FLow
        EE.fx_f(:,end,:) = (repmat(reshape(cell2mat(BC.E(:,1)),1,1,5),size(FF.fx,1),1,1) - EE.fv(:,end,:,2))./GR.dx; % right boundary
    end
    EE.fx_b(:,1,:) = (EE.fv(:,1,:,2) - repmat(reshape(cell2mat(BC.W(:,1)),1,1,5),size(EE.fv,1),1,1))./(GR.dx); % left boundary
    EE.fx = 0.5.*(EE.fx_f + EE.fx_b); % central diff -> extract P_x
    
    %% Calculate Y Derivative
    GG.fy(2:end-1,:,:) = (GG.fv(3:end,:,:) - GG.fv(1:end-2,:,:))./(2*GR.dy);
    GG.fy(end,:,:) = (BC.N{3,1}./BC.N{1,1}.*repmat(reshape(cell2mat(BC.N(:,1)),1,1,5),1,size(GG.fy,2),1) - GG.fv(end-1,:,:))./(2*GR.dy); % multiply v (B/rho) to BC before operating on boundary
    GG_S = repmat(BC.S{3,1},1,1,size(EE.fv,3)).*EE.fv(1,:,:,2); % assume dyBdx = v -> v * EE
    GG_S(:,:,3) = EE.fv(1,:,1,2).*BC.S{3,1}.^2;
    GG.fy(1,:,:) = 0.5.*(diff(GG.fv(1:2,:,:),1,1)./GR.dy + (GG.fv(1,:,:) - GG_S)./(0.5*GR.dy));
    
    EE.fy_f(1:end-1,:,:) = diff(EE.fv(:,:,:,2),1,1)./GR.dy;
    EE.fy_b(2:end,:,:) = EE.fy_f(1:end-1,:,:);
    EE.fy_f(end,:,:) = (repmat(reshape(cell2mat(BC.N(:,1)),1,1,5),1,size(EE.fv,2),1) - EE.fv(end,:,:,2))./GR.dy;
    EE_S = EE.fv(1,:,:,2);
    EE_S(:,:,3) = BC.S{3,1}.* EE.fv(1,:,1,2);
    EE.fy_b(1,:,:) = (EE.fv(1,:,:,2) - EE_S)./(0.5*GR.dy);
    EE.fy = 0.5.*(EE.fy_f + EE.fy_b);
    
    %% Calculate Laplacians 
    
    EE.laplace = (EE.fx_f(:,:,1:4) - EE.fx_b(:,:,1:4))./(GR.dx) + (EE.fy_f(:,:,1:4) - EE.fy_b(:,:,1:4))./(GR.dy);
    
    %% Calculate Next Time Step
    
    DF = (2./(GR.dx^2) + 2./(GR.dy^2)); % Dufort-Frankel coefficient
    pTerms = cat(3, zeros(size(GR.XX)),...
                        EE.fx(:,:,end),...
                        EE.fy(:,:,end),...
                        FF.fx(:,:,end) + GG.fy(:,:,end));
    EE.fv(:,:,1:4,3) = (eps_s.*(EE.laplace + DF.*EE.fv(:,:,1:4,2)) - 0.5*EE.fv(:,:,1:4,1).*(eps_s*DF - 1/dt) - (FF.fx(:,:,1:4) + GG.fy(:,:,1:4) + pTerms))./(0.5/dt + 0.5*eps_s*DF);
    
    %% Calculate Residual
    
    if (size(res,1) == 1) && all(res(end,:) == 0)
        res(1,:) = reshape(max(max(abs(EE.fv(:,:,1:4,3) - EE.fv(:,:,1:4,2)))),1,4,1);
    else
        res(end+1,:) = reshape(max(max(abs(EE.fv(:,:,1:4,3) - EE.fv(:,:,1:4,2)))),1,4,1);
    end

%     figure(1);
% %     hold on;
%     for i = 1:size(res,2)
%         semilogy(1:size(res,1), res(:,i));
%         hold on;
%     end
%     hold off;
%     legend('Density', '\rho u', '\rho v', '\rho \epsilon', 'Location' ,'bestoutside');
    figure(2);
    contourf(GR.XX,GR.YY,EE.fv(:,:,1,3), 50)
    title(['Density (Normalized), M=' num2str(M0)]);
    colorbar('eastoutside');
    axis equal
    drawnow;
    
    if (size(res, 1) > 500) && (mod(size(res, 1), 2000) == 0)
        fprintf('Iteration Ct: %i\n', size(res, 1));
        fprintf('Current Residual: %0.5e\n', max(res(end, :)));
        toc;
        figure(1);semilogy(1:size(res,1), res(:,1));
        hold on;
        semilogy(1:size(res,1), res(:,2));
        semilogy(1:size(res,1), res(:,3));
        hold off;
        legend('Density', '\rho u', '\rho v', 'Location' ,'bestoutside');
%         drawnow;
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

