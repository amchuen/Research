clc;
clear
close all;

%% Initialize Computational Grid and Physical Grid

% Nozzle Function
nozzCoeffs = [2 -2 1]';

% Computational Grid
nZeta = 41; % dZ = 1
nEta = 51; % dE = 1
[ZZ, EE] = meshgrid(1:nZeta, 1:nEta);

% Physical Grid -> boundary conditions
xRange = [0.5, 1];
xS = equalCurveDist(nozzCoeffs, xRange, nZeta);
xW = zeros(nEta, 1);
xN = xS;
xE = ones(nEta, 1);

yN = [0.5.*(1 + (2.*xN(xN <1)-1).^2), ones(size(xN(xN >=1)))];
yS = -yN;
yW = linspace(yS(1), yN(1), nEta);
yE = linspace(yS(end), yN(end), nEta);

UU = cat(3, ZZ, EE); % [X, Y]
for i = 1:size(UU,2)
   UU(:,i,1:2) = cat(3, linspace(xS(i), xN(i), nEta), linspace(yS(i), yN(i), nEta));
end

figure(1);
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
normU = zeros(length(yS),2);
normL = normU;

for i = 1:size(UU,2)
    if i == 1
        normU(i,:) = [-(-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2)),(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))]./sqrt((-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2))^2+(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))^2);
        normL(i,:) = [-(-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2)),(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))]./sqrt((-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2))^2+(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))^2);
    elseif i == size(UU,2)
        normU(i,:) = [-(1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)),(1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2))]./sqrt((1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)).^2 + (1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2)).^2);
        normL(i,:) = [-(1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)),(1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2))]./sqrt((1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)).^2 + (1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2)).^2);
    else
        normU(i,:) = [-(yN(i+1)-yN(i-1)),(xN(i+1)-xN(i-1))]./sqrt((yN(i+1)-yN(i-1))^2+(xN(i+1)-xN(i-1))^2);
        normL(i,:) = [-(yS(i+1)-yS(i-1)),(xS(i+1)-xS(i-1))]./sqrt((yS(i+1)-yS(i-1))^2+(xS(i+1)-xS(i-1))^2);
    end
end

figure(1);hold on;
quiver(xN',yN',normU(:,1), normU(:,2)); hold on;
quiver(xS',yS',normL(:,1), normL(:,2));

%% Fit a Cubic Across y
coeffMat = zeros(4,size(UU,2));
for i = 1:size(UU,2)
    Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
            yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
    bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
    coeffMat(:,i) = Amat\bVec;
    UU(:,i,2) = equalCurveDist(coeffMat(:,i), [yS(i), yN(i)], length(UU(:,i,2)));
    UU(:,i,1) = polyval(coeffMat(:,i), UU(:,i,2));
end

close all;
figure(1);
plot(UU(:,:,1), UU(:,:,2), 'b-'); axis equal;
hold on;
% plot(xN, yN, 'b-', xS, yS, 'b-');

figure();
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-');axis equal;

%% Fit a cubic Across X

% pchipMat = [3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1];
normB = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,1),UU(:,1,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,1),UU(:,1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,1),UU(:,1,2)).^2)];
normMat = zeros(size(UU,1), size(UU,2));
for i = 2:(size(UU,1)-1) % march in y-direction
    
    % Grab Normals at Starting Point
%     y0 = UU(i,1,2);
    y0p = normB(i,2)./normB(i,1);
    
    for ii = 2:size(UU,2) % march in x-direction
        func = @(yVec)  [   yVec(2) + polyval((3:-1:1)'.*coeffMat(1:end-1,ii), yVec(1));...
                            0.5.*(yVec(2)+y0p).*(polyval(coeffMat(:,ii),yVec(1))-UU(i,ii-1,1))+UU(i,ii-1,2)-yVec(1)];
                        
        jacob = @(yVec) [   polyval((2:-1:1)'.*(3:-1:2)'.*coeffMat(1:end-2,ii), yVec(1)), 1;...
                            0.5.*(yVec(2)+y0p).*polyval(coeffMat(1:end-1,ii),yVec(1))-1, 0.5.*(polyval(coeffMat(:,ii),yVec(1))-UU(i,ii-1,1))];
        [yNew, res] = newtonSys(func, jacob, [UU(i,ii,2); -polyval((3:-1:1)'.*coeffMat(1:end-1,ii),UU(i,1,2))], 1e-7, 1e-7, 1e5, 1);
        UU(i,ii,1) = polyval(coeffMat(:,ii),yNew(1)); UU(i,ii,2) = yNew(1);
        y0p = yNew(2);
        tVec = [1, y0p];
        nVec = [-1, polyval((3:-1:1)'.*coeffMat(1:end-1,ii),UU(i,ii,2))];
        normMat(i,ii) = 1-abs((tVec./norm(tVec))*(nVec./norm(nVec))');
    end
    
    figure(1);hold on;
    plot(UU(i,:,1), UU(i,:,2), 'b-');
    drawnow;
    
end

% close all;
figure(); 
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal;
saveas(gcf, 'quadratic_spline_mesh', 'pdf');

figure();contourf(UU(:,:,1), UU(:,:,2), normMat, 50); axis equal;
colorbar;
saveas(gcf, 'quadratic_spline_mesh_orthogCheck', 'pdf');


% figure(1);
% saveas(gcf, 'cubic_mesh_higherRes', 'pdf');