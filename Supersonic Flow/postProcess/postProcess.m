function OUT = postProcess(GR, BC, OUT, dirName)

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

% figure();
% semilogy(1:length(OUT.res), OUT.res(:,1), 1:length(OUT.res), OUT.res(:,2));
% legend('\rho', '\phi','Location', 'Best');
figure(1);
for iii = 1:size(OUT.res,2)
    semilogy(1:length(OUT.res), OUT.res(:,iii));
    hold on;
end
legend(BC.N.varName, 'Location', 'Best');
title('Residual');
xlabel('Iteration');
saveas(gcf, [dirName, 'residual'],'pdf');
saveas(gcf, [dirName, 'residual']);

% Plot Solved Variables
for i = 1:size(OUT.Uvals,3)
    if strcmp(BC.N.varType{i}, 'phi')
        figure();
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end));axis equal;
        title('Flow Potential, \phi');
        xlabel('X');
        ylabel('Y');
        colorbar;
        saveas(gcf, [dirName, 'potential'],'pdf');
        saveas(gcf, [dirName, 'potential']);
        
%         plotCPpotential(GR,BC, OUT.Uvals(:,:,:,end));
    elseif strcmp(BC.N.varType{i}, 's')    
        figure();
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end),50);axis equal;
        title(BC.N.varName{i});
        xlabel('X');
        ylabel('Y');
        colorbar;
        if strcmp(BC.N.varName{i}, '\rho')
            saveas(gcf, [dirName, 'density'],'pdf');
            saveas(gcf, [dirName, 'density']);
        else
            saveas(gcf, [dirName, 'energy'],'pdf');
            saveas(gcf, [dirName, 'energy']);
        end
    
    elseif strcmp(BC.N.varType{i}, 'v1') || strcmp(BC.N.varType{i}, 'v2')
        
        figure();
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end)./OUT.Uvals(:,:,1,end),50);axis equal;
        xlabel('X');
        ylabel('Y');
        colorbar;
        if strcmp(BC.N.varType{i}, 'v2')
            title('V2');
        elseif strcmp(BC.N.varType{i}, 'v1')
            title('V1');
        end
        saveas(gcf, [dirName, BC.N.varType{i}],'pdf');
        saveas(gcf, [dirName, BC.N.varType{i}]);
        
    end
    
end

if any(strcmp(BC.N.varType,'phi'))
    plotCPpotential(OUT.Uvals(:,:,:,end), -10, 10);
elseif any(strcmpi(BC.N.varType, 'v1')) || any(strcmpi(BC.N.varType, 'v2'))
    plotCPeuler(OUT.Uvals(:,:,:,end), -10, 10);
end
    


%end

function plotCPpotential(FF, xstart, xend)
%         FF = OUT.Uvals(:,:,:,end);
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
%             xstart = -4;
%             xend = 4;
            plot(xvals(xvals >=xstart & xvals <= xend), cpX(xvals >=xstart & xvals <= xend));
        else
            xvals = [fliplr(GR.XX(:,end)'),fliplr(GR.XX(1,:))];
            cpX = [fliplr(cpVals(:,end)'),fliplr(cpVals(1,:))];
%             xstart = xvals(1);
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
%         xx = GR.XX(1,2:end-1);
        xx = GR.XX(1,:);
        plot(xx(xx>-2 & xx < 2), cpVals(1,xx>-2 & xx < 2));
        set(gca,'YDir', 'Reverse');
        title('Surface Pressure Coefficient');
        saveas(gcf, [dirName 'cp_surf'],'pdf');
        saveas(gcf, [dirName 'cp_surf']);
    end

end

function plotCPeuler(FF,xstart, xend)
    DIR = 'N';
    indTan = reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
    indNorm = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
    indRho = reshape(strcmp(BC.(DIR).varName, '\rho'),1,1,size(FF,3));
    cpVals = 1 - (FF(:,:,indTan)./FF(:,:,indRho)).^2 - (FF(:,:,indNorm)./FF(:,:,indRho)).^2;
    figure();contourf(GR.XX, GR.YY, cpVals, 50);
    axis equal;
    title('Coefficient of Pressure');
    colorbar;
    saveas(gcf, [dirName 'cp_contour'],'pdf');
    saveas(gcf, [dirName 'cp_contour']);

    figure();
    if GR.isPolar
        if all(GR.YY(:,end) == GR.YY(:,1))
            xvals = [fliplr(GR.XX(:,end)'),fliplr(GR.XX(1,:)), GR.XX(:,1)'];
            cpX = [fliplr(cpVals(:,end)'),fliplr(cpVals(1,:)), cpVals(:,1)'];
%             xstart = -4;
%             xend = 4;
            plot(xvals(xvals >=xstart & xvals <= xend), cpX(xvals >=xstart & xvals <= xend));
        else
            xvals = [fliplr(GR.XX(:,end)'),fliplr(GR.XX(1,:))];
            cpX = [fliplr(cpVals(:,end)'),fliplr(cpVals(1,:))];
%             xstart = xvals(1);
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
        xx = GR.XX(1,:);
%         cpX = cpVals(1,:);
        plot(xx(xx>-2 & xx < 2), cpVals(1,xx>-2 & xx < 2));
        set(gca,'YDir', 'Reverse');
        title('Surface Pressure Coefficient');
        saveas(gcf, [dirName 'cp_surf'],'pdf');
        saveas(gcf, [dirName 'cp_surf']);
    end
    
    PP = FF(:,:,1).^(FL.gam)./(FL.gam .* FL.M0.^2);
    figure();contourf(GR.XX, GR.YY, PP, 50);
    title('Pressure');
    saveas(gcf, [dirName 'pressure'], 'pdf');
    saveas(gcf, [dirName 'pressure']);
    
    figure();plot(GR.XX(1,:), PP(1,:));
    title('Surface Pressure');
    set(gca,'YDir', 'Reverse');
    saveas(gcf, [dirName 'pressure_surf'], 'pdf');
    saveas(gcf, [dirName 'pressure_surf']);

end

end