function OUT = postProcessHalfNozzle(GR, BC, FL, OUT, dirName, eqnFunc)

if ~exist(dirName, 'dir')
    mkdir(dirName);
end

markerArr = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};


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

%% Plot Solved Variables
for i = 1:size(OUT.Uvals,3)
    if strcmp(BC.N.varType{i}, 'phi')
        figure();
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end));%axis equal;
        title('Flow Potential, \phi');
        xlabel('X');
        ylabel('Y');
        colorbar;
        saveas(gcf, [dirName, 'potential'],'pdf');
        saveas(gcf, [dirName, 'potential']);
        
%         plotCPpotential(GR,BC, OUT.Uvals(:,:,:,end));
    elseif strcmp(BC.N.varType{i}, 's')    
        figure();
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end),50);%axis equal;
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
        contourf(GR.XX, GR.YY, OUT.Uvals(:,:,i,end)./OUT.Uvals(:,:,1,end),50);%axis equal;
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

%% Plot Pressure Results

[~, ~, PP] = eqnFunc(GR, FL, BC, OUT.Uvals(:,:,:,end));
figure();
contourf(GR.XX, GR.YY, PP,50);%axis equal;
xlabel('X');
ylabel('Y');
colorbar;
title('Pressure');
saveas(gcf, [dirName, 'pressure']);

%% Plot Centerline Results

figure();
for i = 1:size(OUT.Uvals,3)
    if i == 1
        plot(GR.XX(1,:), OUT.Uvals(end,:,i,end), strcat('-',markerArr{i})); hold on;
    else
        plot(GR.XX(1,:), OUT.Uvals(end,:,i,end)./OUT.Uvals(end,:,1,end), strcat('-',markerArr{i}));
    end
end
plot(GR.XX(1,:), PP(end,:), strcat('-',markerArr{i+1}));
% plot(x_vals, (UU(2,:,3)./UU(1,:,3))./sqrt(gam.*PP./(UU(1,:,3)./g_x)), '--');
% plot(x_vals, 1./sqrt((H0./(UU(2,:,3)./UU(1,:,3)).^2 - 0.5)*(gam-1)), '--');
% title(['Exit Mach Number:' num2str(M_e)]);
legend('\rho', 'u', 'v','E', 'P', 'Location', 'bestoutside');
saveas(gcf, [dirName, 'centerline']);

figure();
plot(GR.XX(1,:), 1./sqrt((3./(OUT.Uvals(end,:,2,3)./OUT.Uvals(end,:,1,3)).^2 - 0.5)*(FL.gam-1)), '--');
title('Mach Number');
xlabel('X');
ylabel('M(x)');
saveas(gcf, [dirName, 'mach_num']);
saveas(gcf, [dirName, 'mach_num'], 'pdf');

end