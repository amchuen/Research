clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
GR.dx = 0.05;
GR.dy = 0.08;

% Field Axis Values
y_max = 7;
x_max = 3 + 20*GR.dx;
x_min = -39*GR.dx; %(-19*dx);
GR.x_vals = x_min:GR.dx:x_max;
GR.y_vals = 0:GR.dy:y_max;
[GR.XX, GR.YY] = meshgrid(GR.x_vals, GR.y_vals);

%% FL - fluid parameters
FL.gam = 1.4; % heat 
% FL.M0 = 1.4;
% FL.M0 = 0.98;

M_inf = [0.95, 0.98, 1.05, 1.1];
alphas = [13, 13, 50, 50];

%% CT - simulation control, including tolerances, viscous factor gain, etc.
CT.v_coeff = 1.5;
CT.tol = 1e-5;
CT.alpha = 13;
CT.dt = 0.05;
CT.t_rho = 0;

%% Boundary Conditions

% % Body Values - Parabolic
tau = 0.05;
x_vals = GR.x_vals;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2.*GR.dx);
end

BC.Vy_II = zeros(size(GR.YY(end,:)));
BC.Vx_II = ones(size(GR.XX(end,:)));
BC.Vx_I = ones(size(GR.XX(:,1)));
BC.PHI_II = GR.XX(end,:);
BC.PHI_I = GR.XX(:,1);
BC.dyBdx = dyBdx;

%% Sim Param Analysis
figure();hold on;

for i = 1:length(M_inf)
    K_sim = (M_inf(i)^2 - 1)/(((FL(1).gam + 1)*M_inf(i)^2 * tau )^(2/3));
    % Step 2 - get M_inf vs tau
%     if M_inf(i) < 1
    del_M = (-0.2*M_inf(i)):0.01:(0.5*M_inf(i));
%     else
%         del_M = (0.0):0.01:(M_inf(i));
%     end

    del_tau = (((M_inf(i)+del_M).^2 - 1).^(1.5))./(K_sim^(1.5) * (FL(1).gam + 1) * (M_inf(i) + del_M).^2) - tau;

    plot(del_M, del_tau);
    
end

legend('M = 0.95', 'M = 0.98', 'M = 1.05', 'M = 1.1', 'Location', 'Best');
xlabel('\Delta M');
ylabel('\Delta \tau');
title('Change in Thickness vs Change in Mach Number of Analogous Fluid');
saveas(gcf, 'sim_param.png');

%% RUN CALCULATION
for i = 1:length(M_inf)
    FL(1).M0 = M_inf(i);
    CT.alpha = alphas(i);
    
    [OUT] = tsf_cart(GR, FL(1), BC, CT);

    %% POST PROCESS

    folderName = ['M_' num2str(FL(1).M0) '_air'];
    geomName = 'sim_param';

    if ~exist([pwd '\' geomName '\' folderName], 'dir')
        mkdir([pwd '\' geomName '\' folderName]);
    end

    close all;
    figure(1);semilogy(1:length(OUT.res), OUT.res);
    title('Residual Plot');
    xlabel('# of iterations');
    ylabel('Residual (Error)');

    saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot']);

    % plot density
    figure();contourf(GR.XX,GR.YY,OUT.rho_ij, 50)
    title('Density (Normalized)');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\density.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\density']);

    figure(); % cp plots
    contourf(GR.XX, GR.YY, 1-(OUT.U_n.^2), 50); %./((RR.*cos(TT)).^2)
    title('Pressure Coefficient Contours');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour']);

    % Surface CP
    figure();
    plot(GR.XX(1,:), 1 - OUT.U_n(1,:).^2);
    xlabel('X');
    %     ylabel('\phi_{\theta}');
    ylabel('C_p');
    title('C_p on surface of Object');
    set(gca, 'Ydir', 'reverse');
    % axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf']);

    figure(); % field potential
    contourf(GR.XX, GR.YY, OUT.PHI(:,:,2), 50);
    title('Field Potential, \Phi');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot']);

    figure(); % theta-dir velocity plots
    contourf(GR.XX, GR.YY, OUT.U_n, 50); %./((RR.*cos(TT)).^2)
    title('\Phi_X velocity');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta']);

    figure();
    contourf(GR.XX, GR.YY, OUT.M_ij, 50);
    title('Mach Number');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\mach.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\mach']);

    %% Save Results
    save([pwd '\' geomName '\' folderName '\results.mat']);


    %% CALCULATE SIMILARITY

    % Step 1 - Get Similarity Param
    K_sim = (FL(1).M0^2 - 1)/(((FL(1).gam + 1)*FL(1).M0^2 * tau )^(2/3));

    % Step 2 - get M_inf vs tau
    del_M = 0.0:0.01:(FL(1).M0);

    del_tau = (((FL(1).M0+del_M).^2 - 1).^(1.5))./(K_sim^(1.5) * (FL(1).gam + 1) * (FL(1).M0 + del_M).^2) - tau;

    figure();plot(del_M, del_tau);

    % Step 3 - Calculate Similarity for Water
    % Iterate to get the new Mach Number
    FL(2) = FL(1);
    FL(2).gam = 2.0; % water

    M1 = 0.98;
    M2 = sqrt(1 + K_sim*((FL(2).gam +1).*M1^2 * tau)^(2/3));

    res = abs((M2 - M1)/M1);

    while res > 1e-5
       M1 = M2;
       M2 = sqrt(1 + K_sim*(((FL(2).gam +1).*M1^2 * tau)^(2/3)));

       res = abs((M2 - M1)/M1);
    end

    fprintf('Equivalent Mach Number for Water: %0.5f\n', M2);

    FL(2).M0 = M2;

    K_wat = (FL(2).M0^2 - 1)/(((FL(2).gam + 1)*FL(2).M0^2 * tau )^(2/3));

    %% RUN WATER RESULTS

    [OUT_w] = tsf_cart(GR, FL(2), BC, CT);

    folderName = ['M_' num2str(FL(2).M0) '_water'];
    geomName = 'sim_param';

    if ~exist([pwd '\' geomName '\' folderName], 'dir')
        mkdir([pwd '\' geomName '\' folderName]);
    end

    close all;
    figure(1);semilogy(1:length(OUT_w.res), OUT_w.res);
    title('Residual Plot');
    xlabel('# of iterations');
    ylabel('Residual (Error)');

    saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\residual_plot']);

    % plot density
    figure();contourf(GR.XX,GR.YY,OUT_w.rho_ij, 50)
    title('Density (Normalized)');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\density.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\density']);

    figure(); % cp plots
    contourf(GR.XX, GR.YY, 1-(OUT_w.U_n.^2), 50); %./((RR.*cos(TT)).^2)
    title('Pressure Coefficient Contours');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_contour']);

    % Surface CP
    figure();
    plot(GR.XX(1,:), 1 - OUT_w.U_n(1,:).^2);
    xlabel('X');
    %     ylabel('\phi_{\theta}');
    ylabel('C_p');
    title('C_p on surface of Object');
    set(gca, 'Ydir', 'reverse');
    % axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\cp_surf']);

    figure(); % field potential
    contourf(GR.XX, GR.YY, OUT_w.PHI(:,:,2), 50);
    title('Field Potential, \Phi');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_pot']);

    figure(); % theta-dir velocity plots
    contourf(GR.XX, GR.YY, OUT_w.U_n, 50); %./((RR.*cos(TT)).^2)
    title('\Phi_X velocity');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\phi_theta']);

    figure();
    contourf(GR.XX, GR.YY, OUT_w.M_ij, 50);
    title('Mach Number');
    colorbar('eastoutside');
    axis equal
    saveas(gcf, [pwd '\' geomName '\' folderName '\mach.png']);
    saveas(gcf, [pwd '\' geomName '\' folderName '\mach']);

    %% Save Results
    save([pwd '\' geomName '\' folderName '\results.mat']);
    
    clearvars -except i M_inf alphas CT FL BC GR tau
end
