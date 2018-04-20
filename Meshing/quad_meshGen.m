clc;
clear
close all;

%% Initialize Computational Grid and Physical Grid

% Nozzle Function
nozzCoeffs = [2 -2 1]';

% Computational Grid
nZeta = 11; % dZ = 1
nEta = 21; % dE = 1
[ZZ, EE] = meshgrid(1:nZeta, 1:nEta);

% Physical Grid -> boundary conditions
xRange = [0.5, 1];
xS = equalCurveDist(nozzCoeffs, xRange, 2*(nZeta-1)+1);
xN = xS;
xW = linspace(xS(1), xN(1), 2*(nEta-1)+1);%, 1);
xE = ones(2*(nEta-1)+1, 1);

yN = [0.5.*(1 + (2.*xN(xN <1)-1).^2), ones(size(xN(xN >=1)))];
yS = -yN;
yW = linspace(yS(1), yN(1), 2*(nEta-1)+1);
yE = linspace(yS(end), yN(end), 2*(nEta-1)+1);

% Generate Primary Grid
UU_1 = cat(3, ZZ, EE); % [X, Y]
normU = zeros(length(yS),2);
normL = normU;
coeffMat_1 = zeros(4,size(UU_1,2));
for i = 1:size(xN,2)
%     ind = 2*i-1;
%     UU(:,i,1:2) = cat(3, linspace(xS(ind), xN(ind), nEta), linspace(yS(ind), yN(ind), nEta));
%     if i == 1
%         normU(i,:) = [-(-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2)),(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))]./sqrt((-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2))^2+(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))^2);
%         normL(i,:) = [-(-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2)),(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))]./sqrt((-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2))^2+(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))^2);
%     elseif i == size(xN,2) %size(UU,2)
%         normU(i,:) = [-(1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)),(1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2))]./sqrt((1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)).^2 + (1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2)).^2);
%         normL(i,:) = [-(1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)),(1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2))]./sqrt((1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)).^2 + (1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2)).^2);
%     else
%         normU(i,:) = [-(yN(i+1)-yN(i-1)),(xN(i+1)-xN(i-1))]./sqrt((yN(i+1)-yN(i-1))^2+(xN(i+1)-xN(i-1))^2);
%         normL(i,:) = [-(yS(i+1)-yS(i-1)),(xS(i+1)-xS(i-1))]./sqrt((yS(i+1)-yS(i-1))^2+(xS(i+1)-xS(i-1))^2);
%     end
    normU(i,:) = [-2*(2*xN(i)-1), 1]./sqrt((-2*(2*xN(i)-1)).^2 + 1);
    normL(i,:) = [2*(2*xN(i)-1), 1]./sqrt((2*(2*xN(i)-1)).^2 + 1);
    
    if mod(i,2) ~= 0 % update primary grid
        ind = 0.5*(i + 1);
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_1(:,ind) = Amat\bVec;
        UU_1(:,ind,2) = equalCurveDist(coeffMat_1(:,ind), [yS(i), yN(i)], length(UU_1(:,ind,2)));
        UU_1(:,ind,1) = polyval(coeffMat_1(:,ind), UU_1(:,ind,2));
    end

end

% build secondary grid
UU_2 = zeros(size(UU_1));
UU_2(end+1, end+1, :) = 0;
UU_2(:, 1,1) = xW([1,2:2:length(xW)-1,length(xW)]);
UU_2(:, 1,2) = yW([1,2:2:length(xW)-1,length(xW)]);
UU_2(end, :,1) = xN([1,2:2:length(xN)-1,length(xN)]);
UU_2(end, :,2) = yN([1,2:2:length(yN)-1,length(yN)]);
UU_2(1, :,1) = xS([1,2:2:length(xS)-1,length(xS)]);
UU_2(1, :,2) = yS([1,2:2:length(yS)-1,length(yS)]);
coeffMat_2 = zeros(4,size(UU_2,2));
for i = [1,2:2:length(xN)-1,length(xN)] %2:2:(size(xN,2)-1)
    if i == 1
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,i) = Amat\bVec;
%         UU_2(:,i,2) = equalCurveDist(coeffMat_2(:,i), [yS(i), yN(i)], length(UU_2(:,i,2)));
%         UU_2(:,i,1) = polyval(coeffMat_2(:,i), UU_2(:,i,2));
    elseif i == length(xN)
        y_S = 0.5*(UU_1(1,end,2) + UU_1(2,end,2));
        y_N = 0.5*(UU_1(end,end,2) + UU_1(end-1,end,2));
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,end) = Amat\bVec;
        UU_2(2:end-1,end,2) = equalCurveDist(coeffMat_2(:,end), [y_S, y_N], length(UU_2(2:end-1,end,2)));
        UU_2(2:end-1,end,1) = polyval(coeffMat_2(:,end), UU_2(2:end-1,end,2));
    else
        ind = 0.5*(i);
        y_S = 0.25*(UU_1(1,ind,2) + UU_1(1,ind+1,2) + UU_1(2,ind,2) + UU_1(2,ind+1,2));
        y_N = 0.25*(UU_1(end,ind,2) + UU_1(end,ind+1,2) + UU_1(end-1,ind,2) + UU_1(end-1,ind+1,2));
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,ind+1) = Amat\bVec;
        UU_2(2:end-1,ind+1,2) = equalCurveDist(coeffMat_2(:,ind+1), [y_S, y_N], length(UU_2(2:end-1,ind+1,2)));
        UU_2(2:end-1,ind+1,1) = polyval(coeffMat_2(:,ind+1), UU_2(2:end-1,ind+1,2));
    end
    
end

figure(1);hold on;
plot(UU_1(:,:,1), UU_1(:,:,2), 'b-', UU_1(:,:,1)', UU_1(:,:,2)', 'b-'); axis equal; 
hold on;
plot(UU_2(:,:,1), UU_2(:,:,2), 'r-', UU_2(:,:,1)', UU_2(:,:,2)', 'r-');
quiver(xN(1:2:end)',yN(1:2:end)',normU(1:2:end,1), normU(1:2:end,2)); hold on;
quiver(xS(1:2:end)',yS(1:2:end)',normL(1:2:end,1), normL(1:2:end,2));
plot(xN, yN, 'o');
plot(xS, yS, 'o');

%% Fit a quadratic Across X

% Refit Primary Grid
normB_1 = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2)).^2), polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2)).^2)];
normMat_1 = zeros(size(UU_1,1), size(UU_1,2));
for i = 2:(size(UU_1,1)-1) % march in y-direction
    
    % Grab Normals at Starting Point
    y0p = normB_1(i,2)./normB_1(i,1);
    
    for ii = 2:size(UU_1,2) % march in x-direction
        func = @(yVec)  [   yVec(2) + polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii), yVec(1));...
                            0.5.*(yVec(2)+y0p).*(polyval(coeffMat_1(:,ii),yVec(1))-UU_1(i,ii-1,1))+UU_1(i,ii-1,2)-yVec(1)];
                        
        jacob = @(yVec) [   polyval((2:-1:1)'.*(3:-1:2)'.*coeffMat_1(1:end-2,ii), yVec(1)), 1;...
                            0.5.*(yVec(2)+y0p).*polyval(coeffMat_1(1:end-1,ii),yVec(1))-1, 0.5.*(polyval(coeffMat_1(:,ii),yVec(1))-UU_1(i,ii-1,1))];
        [yNew, res] = newtonSys(func, jacob, [UU_1(i,ii,2); -polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii),UU_1(i,1,2))], 1e-7, 1e-7, 1e5, 1);
        UU_1(i,ii,1) = polyval(coeffMat_1(:,ii),yNew(1)); UU_1(i,ii,2) = yNew(1);
        y0p = yNew(2);
        tVec = [1, y0p];
        nVec = [-1, polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii),UU_1(i,ii,2))];
        normMat_1(i,ii) = 1-abs((tVec./norm(tVec))*(nVec./norm(nVec))');
    end
    
end

% Refit Secondary Grid
normB_2 = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2)).^2), polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2)).^2)];
normMat_2 = zeros(size(UU_2,1), size(UU_2,2));
for i = 2:(size(UU_2,1)-1) % march in y-direction
    
    % Grab Normals at Starting Point
    y0p = normB_1(i,2)./normB_1(i,1);
    
    for ii = 2:size(UU_2,2) % march in x-direction
        func = @(yVec)  [   yVec(2) + polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii), yVec(1));...
                            0.5.*(yVec(2)+y0p).*(polyval(coeffMat_2(:,ii),yVec(1))-UU_2(i,ii-1,1))+UU_2(i,ii-1,2)-yVec(1)];
                        
        jacob = @(yVec) [   polyval((2:-1:1)'.*(3:-1:2)'.*coeffMat_2(1:end-2,ii), yVec(1)), 1;...
                            0.5.*(yVec(2)+y0p).*polyval(coeffMat_2(1:end-1,ii),yVec(1))-1, 0.5.*(polyval(coeffMat_2(:,ii),yVec(1))-UU_2(i,ii-1,1))];
        [yNew, res] = newtonSys(func, jacob, [UU_2(i,ii,2); -polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii),UU_2(i,1,2))], 1e-7, 1e-7, 1e5, 1);
        UU_2(i,ii,1) = polyval(coeffMat_2(:,ii),yNew(1)); UU_2(i,ii,2) = yNew(1);
        y0p = yNew(2);
        tVec = [1, y0p];
        nVec = [-1, polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii),UU_2(i,ii,2))];
        normMat_2(i,ii) = 1-abs((tVec./norm(tVec))*(nVec./norm(nVec))');
    end
    
end

% close all;
figure(10); 
plot(UU_1(:,:,1), UU_1(:,:,2), 'b-', UU_1(:,:,1)', UU_1(:,:,2)', 'b-'); axis equal;hold on;
plot(UU_2(:,:,1), UU_2(:,:,2), 'r--', UU_2(:,:,1)', UU_2(:,:,2)', 'r--'); axis equal;
saveas(gcf, 'quadratic_spline_mesh', 'pdf');

% figure();contourf(UU_1(:,:,1), UU_1(:,:,2), normMat_1, 50); axis equal;
% colorbar;
% saveas(gcf, 'quadratic_spline_mesh_orthogCheck', 'pdf');


% figure(1);
% saveas(gcf, 'cubic_mesh_higherRes', 'pdf');

% %% Find Cell Volumes
% deltaX_NS = UU_1(:,2:end,1:2) - UU_1(:,1:end-1,1:2);
% xTest = cat(3,1,-1).*deltaX_NS(:,:,2:-1:1);%[deltaY, -deltaX];%./sqrt(deltaY.^2 + deltaX.^2);
% midX_NS = 0.5.*(UU_1(:,2:end,:) + UU_1(:,1:end-1,:));
% xx1 = -xTest(1,1,:) + midX_NS(1,1,:);
% xx2 = xTest(2,1,:) + midX_NS(2,1,:);
% figure(10); hold on;
% plot([midX_NS(1,1,1), xx1(1)], [midX_NS(1,1,2), xx1(2)], 'ro-');
% plot([midX_NS(2,1,1), xx2(1)], [midX_NS(2,1,2), xx2(2)], 'ko-');
% 
% deltaX_EW = UU_1(2:end,:,:) - UU_1(1:end-1,:,:);
% xTest_EW = cat(3,1,-1).*deltaX_EW(:,:,2:-1:1);
% midX_EW = 0.5.*(UU_1(2:end,:,:) + UU_1(1:end-1,:,:));
% xx3 = xTest_EW(1,1,:) + midX_EW(1,1,:);
% xx4 = -xTest_EW(1,2,:) + midX_EW(1,2,:);
% 
% plot([midX_EW(1,1,1), xx3(1)], [midX_EW(1,1,2), xx3(2)], 'bo-');
% plot([midX_EW(1,2,1), xx4(1)], [midX_EW(1,2,2), xx4(2)], 'o-');