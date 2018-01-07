function post_process(OUT, GR, geomName, folderName)

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
% string = sprintf('save([pwd ''\\'' %s ''\\'' %s ''\\results.mat'']);', geomName, folderName);
string = 'save([pwd ''\'' geomName ''\'' folderName ''\results.mat'']);';

evalin('base', string);

end