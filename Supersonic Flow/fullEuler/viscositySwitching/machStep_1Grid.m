clc;
close all;
clear;

%% CT - simulation control, including tolerances, viscous factor gain, etc.
eps_s = 1e-4; % spatial diffusion term
% eps_t = 0.005; % time diffusion term
tol = 1e-5;
dt = 2.16e-3;
iter_min = 300;
CFL_on = 1;

case_name = 'machStep';
    
%% FL - fluid parameters
gam = 1.4; % heat 
M0 = 3.0;

%% GR - grid information, such as the meshfield, grid spacing (dr, dT, etc.)

% Define Grid
dx = 0.01;
dy = 0.01;

% Field Axis Values - body fitted grid
% Field Axis Values - body fitted grid
% x_range=[   -4-39*dx,...
%             7 + 20*dx];
x_range=[   0.5*dx,...
            3];
y_range=[   0.5*dy,...
            1-0.5*dy];
x_wall = 0.6;
y_wall = 0.2;

GR = struct('Y_vals',[], 'X_vals',[], 'XX',[], 'YY',[]);

GR.y_vals = y_range(1):dy:y_range(2);
GR.x_vals = x_range(1):dx:x_range(2);
GR.dy = dy; GR.dx = dx;

[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

% Show Grid
figure();
plot(GR.XX,GR.YY, 'b-', GR.XX', GR.YY', 'b-');
axis equal
hold on;
rectangle('Position',[x_wall 0 3-x_wall y_wall], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
stepInd = {GR.y_vals < y_wall, GR.x_vals > x_wall};
wallIndX = (length(GR.x_vals) - sum(stepInd{2})+1);
        
%% INITIALIZE

EE = struct('fv',[],'GR', []);
FF = EE; 
GG = EE; 

% Matrix Dimensions
% 1:Y, 2:X, 3:vec, 4:t

EE.fv = repmat(cat(3,   ones(size(GR.XX)),... % rho
                        ones(size(GR.XX)),... % rho * u
                        zeros(size(GR.XX)),...% rho * v
                        ones(size(GR.XX))./(gam*(gam-1)*M0^2) + 0.5,... % rho * E
                        1/(gam * M0^2).*ones(size(GR.XX))),... % P ... eps_mach error here?
                        1, 1, 1, 3); 
EE.fv(stepInd{1}, stepInd{2},:,:) = NaN;

% EE.fv(:,:,5,:) = repmat((EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1),1,1,1,3);
EE.fx_f = zeros(size(EE.fv(:,:,:,2)));
EE.fx_b = EE.fx_f;
EE.fy_f = EE.fx_f;
EE.fy_b = EE.fx_f;
EE.eps1f = EE.fx_f;
EE.eps1b = EE.fx_f;
EE.avg1f = EE.fx_f;
EE.avg1b = EE.fx_f;
EE.eps2f = EE.fx_f;
EE.eps2b = EE.fx_f;
EE.avg2f = EE.fx_f;
EE.avg2b = EE.fx_f;

FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
FF.fv(stepInd{1}, stepInd{2},:) = NaN;
FF.fx = zeros(size(FF.fv));

GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
GG.fv(stepInd{1}, stepInd{2},:) = NaN;
GG.fy = zeros(size(GG.fv));

%% Boundary Condition

% Vector E:
BC.W = {    1.,  'D';... rho
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
            0, 'D';... % wall is normal to this direction, i.e. wall is parallel to x-dir
            0,  'N';...
            0,  'N',...
            };

%% Solve

res = [];
tic;
tot_time = 0;

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
    
    % manage wall step BC - reflective wall
    FF.fv(stepInd{1}, wallIndX, :) = -FF.fv(stepInd{1}, wallIndX-1, :);
    FF.fv(stepInd{1}, wallIndX, 2) = -FF.fv(stepInd{1}, wallIndX, 2);
    EE.fv(stepInd{1}, wallIndX, :,2) = EE.fv(stepInd{1}, wallIndX-1, :,2);
    EE.fv(stepInd{1}, wallIndX, 2,2) = -EE.fv(stepInd{1}, wallIndX, 2,2);
    
    FF.fx(:,2:end-1,:) = (FF.fv(:,3:end,:) - FF.fv(:,1:end-2,:))./(2*GR.dx); % central difference
    FF.fx(:,end,:) = repmat(reshape(cell2mat(BC.E(:,1)),1,1,5),size(FF.fx,1),1,1);
    FF.fx(:,1,:) = (FF.fv(:,2,:) - BC.W{2,1}./BC.W{1,1}.*repmat(reshape(cell2mat(BC.W(:,1)),1,1,5),size(FF.fx,1),1,1))./(2*GR.dx);

    EE.fx_f(:,1:end-1,:) = diff(EE.fv(:,:,:,2),1,2)./GR.dx; % forward difference
    EE.avg1f(:,1:end-1,1) = 0.5.*(abs(EE.fv(:,2:end,1,2))+abs(EE.fv(:,1:end-1,1,2))); % averaging denominator?
    EE.fx_b(:,2:end,:) = EE.fx_f(:,1:end-1,:); % backward difference
    EE.avg1b(:,2:end,1) = EE.avg1f(:,1:end-1,1);
    EE.fx_f(:,end,:) = (EE.fv(:,end-1,:,2) - EE.fv(:,end,:,2))./GR.dx; % right boundary, 2nd order Neumann condition
    EE.avg1f(:,end,1) = 0.5.*(abs(EE.fv(:,end-1,1,2))+abs(EE.fv(:,end,1,2)));
    EE.fx_b(:,1,:) = (EE.fv(:,1,:,2) - repmat(reshape(cell2mat(BC.W(:,1)),1,1,5),size(EE.fv,1),1,1))./(GR.dx); % left boundary
    EE.avg1b(:,1,1) = 0.5.*(abs(EE.fv(:,1,1,2)) + abs(BC.W{1,1}));
    EE.fx = 0.5.*(EE.fx_f + EE.fx_b); % central diff -> extract P_x
    EE.eps1f = repmat((EE.fx_f(:,:,1)).^2./EE.avg1f(:,:,1), 1,1,4);
    EE.eps1b = repmat((EE.fx_b(:,:,1)).^2 ./ EE.avg1b(:,:,1),1,1,4);
    
    %% Calculate Y Derivative
    
    % manage wall step BC - reflective wall?
    GG.fv(sum(stepInd{1}), stepInd{2}, :) = -GG.fv(sum(stepInd{1})+1, stepInd{2}, :);
    GG.fv(sum(stepInd{1}), stepInd{2}, 3) = -GG.fv(sum(stepInd{1}), stepInd{2}, 3);
    EE.fv(sum(stepInd{1}), stepInd{2}, :,2) = EE.fv(sum(stepInd{1})+1, stepInd{2}, :,2);
    EE.fv(sum(stepInd{1}), stepInd{2}, 3,2) = -EE.fv(sum(stepInd{1}), stepInd{2}, 3,2);
    
    GG.fy(2:end-1,:,:) = (GG.fv(3:end,:,:) - GG.fv(1:end-2,:,:))./(2*GR.dy);
    GG.fy(end,:,:) = (-GG.fv(end,:,:) - GG.fv(end-1,:,:))./(2*GR.dy); % change sign b/c of reflected boundary
    GG.fy(end,:,3) = (GG.fv(end,:,3) - GG.fv(end-1,:,3))./(2*GR.dy); % scalar vals are not reflected in direction
    GG.fy(1,:,:) = (GG.fv(2,:,:) + GG.fv(1,:,:))./(2*GR.dy);
    GG.fy(1,:,3) = (GG.fv(2,:,3) - GG.fv(1,:,3))./(2*GR.dy); % scalar vals are not reflected in direction
    
    EE.fy_f(1:end-1,:,:) = diff(EE.fv(:,:,:,2),1,1)./GR.dy;
    EE.avg2f(1:end-1,:,:) = 0.5.*(abs(EE.fv(2:end,:,:,2))+abs(EE.fv(1:end-1,:,:,2)));
    EE.fy_b(2:end,:,:) = EE.fy_f(1:end-1,:,:);
    EE.avg2b(2:end,:,:) = EE.avg2f(1:end-1,:,:);
    EE.fy_f(end,:,:) = (EE.fv(end,:,:,2) - EE.fv(end,:,:,2))./GR.dy;
    EE.fy_f(end,:,3) = (-EE.fv(end,:,3,2) - EE.fv(end,:,3,2))./GR.dy;
    EE.avg2f(end,:,:) = 0.5.*(abs(EE.fv(end,:,:,2)) + abs(EE.fv(end,:,:,2)));
    EE.fy_b(1,:,:) = (EE.fv(1,:,:,2) - EE.fv(1,:,:,2))./(GR.dy);
    EE.fy_b(1,:,3) = (EE.fv(1,:,3,2) + EE.fv(1,:,3,2))./(GR.dy); % direction must be reflected across boundary
    EE.avg2b(1,:,:) = 0.5.*(abs(EE.fv(1,:,:,2)) + abs(EE.fv(1,:,:,2)));
    EE.fy = 0.5.*(EE.fy_f + EE.fy_b);
    EE.eps2f = repmat((EE.fy_f(:,:,1)).^2 ./ EE.avg2f(:,:,1), 1,1,4);
    EE.eps2b = repmat((EE.fy_b(:,:,1)).^2 ./ EE.avg2b(:,:,1),1,1,4);
    
    %% Calculate Laplacians 
    
    EE.laplace = (EE.eps1f.*EE.fx_f(:,:,1:4) - EE.eps1b.*EE.fx_b(:,:,1:4))./(GR.dx) + (EE.eps2f.*EE.fy_f(:,:,1:4) - EE.eps2b.*EE.fy_b(:,:,1:4))./(GR.dy);
    
    %% Calculate Next Time Step
    
    DF = ((EE.eps1f+EE.eps1b)./(GR.dx^2) + (EE.eps2f+EE.eps2b)./(GR.dy^2)); % Dufort-Frankel coefficient
    pTerms = cat(3, zeros(size(GR.XX)),...
                        EE.fx(:,:,end),...
                        EE.fy(:,:,end),...
                        FF.fx(:,:,end) + GG.fy(:,:,end));
    EE.fv(:,:,1:4,3) = (eps_s.*(EE.laplace + DF.*EE.fv(:,:,1:4,2)) - 0.5*EE.fv(:,:,1:4,1).*(eps_s*DF - 1/dt) - (FF.fx(:,:,1:4) + GG.fy(:,:,1:4) + pTerms))./(0.5/dt + 0.5*eps_s*DF);
    EE.fv(stepInd{1}, stepInd{2},:,:) = NaN;
    tot_time = tot_time + dt;
    
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
    title(sprintf('Density (Normalized), t=%0.5e, dx=%0.2f', tot_time, dx));
    colorbar('eastoutside');
    axis equal
    drawnow;
    if (abs(4 - tot_time) < dt) && (tot_time > 4)
        saveas(gcf, [case_name '_' num2str(dx) '_' num2str(eps_s) '_density.pdf']);
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
        legend('Density', '\rho u', '\rho v', 'Location' ,'bestoutside');
%         drawnow;
        fprintf('\n');
    end
    
    if (~isreal(EE.fv(not(stepInd{1}),not(stepInd{2}),:,end)) || any(any(any(isnan(EE.fv(not(stepInd{1}),not(stepInd{2}),:,end))))))
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

