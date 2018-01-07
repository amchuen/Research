function postProcess(GR, OUT, dirName)

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

figure();
semilogy(1:length(OUT.res), OUT.res(:,1), 1:length(OUT.res), OUT.res(:,2));
legend('\rho', '\phi','Location', 'Best');
title('Residual');
xlabel('Iteration');
saveas(gcf, [dirName, 'residual'],'pdf');
saveas(gcf, [dirName, 'residual']);

figure();
contourf(GR.XX, GR.YY, OUT.Uvals(:,:,2,end));axis equal;
title('Flow Potential, \phi');
xlabel('X');
ylabel('Y');
colorbar;
saveas(gcf, [dirName, 'potential'],'pdf');
saveas(gcf, [dirName, 'potential']);

figure();
contourf(GR.XX, GR.YY, OUT.Uvals(:,:,1,end),200);axis equal;
title('Density, \rho');
xlabel('X');
ylabel('Y');
colorbar;
saveas(gcf, [dirName, 'density'],'pdf');
saveas(gcf, [dirName, 'density']);

phiX = (diff(OUT.Uvals(:,1:end-1,2,end),1,2) + diff(OUT.Uvals(:,2:end,2,end),1,2))./(2.*GR.dx);
phiYs = (OUT.Uvals(2,2:end-1,2,end) - OUT.Uvals(1,2:end-1,2,end))./GR.dy;
figure();contourf(GR.XX(:,2:end-1), GR.YY(:,2:end-1),phiX, 200);
axis equal;
title('Velocity, U');
colorbar;
saveas(gcf, [dirName 'uVel'],'pdf');
saveas(gcf, [dirName 'uVel']);

figure();
xx = GR.XX(1,2:end-1);
plot(xx(xx>-2 & xx < 2), 1-phiX(1,xx>-2 & xx < 2).^2-phiYs(xx>-2 & xx < 2).^2);
set(gca,'YDir', 'Reverse');
title('Surface Pressure Coefficient');
saveas(gcf, [dirName 'cp_surf'],'pdf');
saveas(gcf, [dirName 'cp_surf']);

end