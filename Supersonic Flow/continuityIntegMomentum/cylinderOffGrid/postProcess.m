function postProcess(GR, BC, OUT, dirName)

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
contourf(GR.XX, GR.YY, OUT.Uvals(:,:,1,end),50);axis equal;
title('Density, \rho');
xlabel('X');
ylabel('Y');
colorbar;
saveas(gcf, [dirName, 'density'],'pdf');
saveas(gcf, [dirName, 'density']);

FF = OUT.Uvals(:,:,:,end);
if GR.isPolar
    phi2_w = diff([bcCalc(GR, BC, FF, 'W', 2), FF(:,:,2)],1,2)./(GR.RR.*GR.dT);
    phi2_e = diff([FF(:,:,2), bcCalc(GR, BC, FF, 'E', 2)],1,2)./(GR.RR.*GR.dT);
    phi1_n = diff([FF(:,:,2); bcCalc(GR, BC, FF, 'N', 2)],1,1)./GR.dR;
    phi1_s = diff([bcCalc(GR, BC, FF, 'S', 2); FF(:,:,2)],1,1)./GR.dR; 
    
    %fvOUT(:,:,1) = (rho_e.*phi2_e - rho_w.*phi2_w)./(GR.dT.*GR.RR) + (GR.RR_N.*rho_n.*phi1_n - GR.RR_S.*rho_s.*phi1_s)./(GR.dR .* GR.RR);
else
    phi2_w = diff([bcCalc(GR, BC, FF, 'W', 2), FF(:,:,2)],1,2)./GR.dx;
    phi2_e = diff([FF(:,:,2), bcCalc(GR, BC, FF, 'E', 2)],1,2)./GR.dx;
    phi1_n = diff([FF(:,:,2); bcCalc(GR, BC, FF, 'N', 2)],1,1)./GR.dy;
    phi1_s = diff([bcCalc(GR, BC, FF, 'S', 2); FF(:,:,2)],1,1)./GR.dy;
    
    %fvOUT(:,:,1) = (rho_e.*phi2_e - rho_w.*phi2_w)./GR.dx + (rho_n.*phi1_n - rho_s.*phi1_s)./GR.dy;
end

phi2 = 0.5.*(phi2_w + phi2_e);
phi1 = 0.5.*(phi1_n + phi1_s);
cpVals = 1-(phi2.^2 + phi1.^2);
figure();contourf(GR.XX, GR.YY, cpVals, 50);
axis equal;
title('Coefficient of Pressure, U');
colorbar;
saveas(gcf, [dirName 'cp_contour'],'pdf');
saveas(gcf, [dirName 'cp_contour']);

figure();
if GR.isPolar
    if all(GR.YY(:,end) == GR.YY(:,1))
        xvals = [fliplr(GR.XX(:,end)'),fliplr(GR.XX(1,:)), GR.XX(:,1)'];
        cpX = [fliplr(cpVals(:,end)'),fliplr(cpVals(1,:)), cpVals(:,1)'];
        xstart = -4;
        xend = 4;
        plot(xvals(xvals >=xstart & xvals <= xend), cpX(xvals >=xstart & xvals <= xend));
    else
        xvals = [fliplr(GR.XX(:,end)'),fliplr(GR.XX(1,:))];
        cpX = [fliplr(cpVals(:,end)'),fliplr(cpVals(1,:))];
        xstart = xvals(1);
        plot(xvals(xvals >=xstart), cpX(xvals >=xstart));
    end
    set(gca,'YDir', 'Reverse');
    title('Surface Pressure Coefficient - X axis');
    saveas(gcf, [dirName 'cp_surfX'],'pdf');
    saveas(gcf, [dirName 'cp_surfX']);
    
    figure();
    plot(GR.TT(1,:), cpVals(1,:));
    set(gca,'XDir', 'Reverse');
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'}) 
    set(gca,'YDir', 'Reverse');
    title('Surface Pressure Coefficient');
    saveas(gcf, [dirName 'cp_surf'],'pdf');
    saveas(gcf, [dirName 'cp_surf']);
else
    xx = GR.XX(1,2:end-1);
    plot(xx(xx>-2 & xx < 2), 1-phiX(1,xx>-2 & xx < 2).^2-phiYs(xx>-2 & xx < 2).^2);
    set(gca,'YDir', 'Reverse');
    title('Surface Pressure Coefficient');
    saveas(gcf, [dirName 'cp_surf'],'pdf');
    saveas(gcf, [dirName 'cp_surf']);
end


end