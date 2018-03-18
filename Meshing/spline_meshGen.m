clc;
clear
close all;

%% Control Parameters
tol = 1e-5;
omega = 1; % change later when norms are easier to find?
cBC = 1/20;

%% Initialize Computational Grid and Physical Grid

% Nozzle Function
nozzCoeffs = [2 -2 1]';

% Computational Grid
nZeta = 25; % dZ = 1
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
plot(xN, yN, 'b-', xS, yS, 'b-');

%% Fit a cubic Across X between each grid line
% pchipMat = [3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1];
% for i = 2:size(UU,2) % march in x-direction
%     
%     if i == 2
%         normB = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2)).^2)];
%     else
%         normB = normF;
%     end
%     
%     normF = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2)).^2)];
%     deltaX = UU(:,i,1) - UU(:,i-1,1);
%     for ii = 1:size(UU,1) %2:(size(UU,1)-1) % march in y-direction
%         coeffs = pchipMat\[deltaX(ii)*normF(ii,2)/normF(ii,1); deltaX(ii)*normB(ii,2)/normB(ii,1); UU(ii,i,2); UU(ii,i-1,2)];
%         xFit = linspace(UU(ii,i-1,1),UU(ii,i,1),11);
%         pFit = polyval(coeffs, (xFit-UU(ii,i-1,1))./deltaX(ii));
%         figure(1);hold on;
%         plot(xFit, pFit, 'b-');
%         drawnow;
%     end
%     
% end

pchipMat = [2 1 0; 0 1 0; 0 0 1];
for i = 2:size(UU,2) % march in x-direction
    
    
    normB = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i-1),UU(:,i-1,2)).^2)];    
    normF = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),UU(:,i,2)).^2)];
    deltaX = UU(:,i,1) - UU(:,i-1,1);
    if i == 6
       fprintf('stop here\n'); 
    end
    for ii = 2:(size(UU,1)-1) %2:(size(UU,1)-1) % march in y-direction
        xOld = UU(ii,i,1);
        coeffs = pchipMat\[deltaX(ii)*normF(ii,2)/normF(ii,1); deltaX(ii)*normB(ii,2)/normB(ii,1); UU(ii,i-1,2)];
        xNew = polyval(coeffMat(:,i),sum(coeffs));
        normF(ii,:) = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs)).^2)];
        res = abs(xNew - xOld);
        
        while res(end) > 1e-5
            xOld = xNew;
            deltaX(ii) = xOld - UU(ii,i-1,1);
            coeffs = pchipMat\[deltaX(ii)*normF(ii,2)/normF(ii,1); deltaX(ii)*normB(ii,2)/normB(ii,1); UU(ii,i-1,2)];
            xNew = polyval(coeffMat(:,i),sum(coeffs));
            normF(ii,:) = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,i),sum(coeffs)).^2)];
            res(end+1) = abs(xNew - xOld);
        end
        
        deltaX(ii) = xNew - UU(ii,i-1,1);
        coeffs = pchipMat\[deltaX(ii)*normF(ii,2)/normF(ii,1); deltaX(ii)*normB(ii,2)/normB(ii,1); UU(ii,i-1,2)];
        UU(ii,i,1) = xNew;
        UU(ii,i,2) = sum(coeffs);
        xFit = linspace(UU(ii,i-1,1),UU(ii,i,1),11);
        pFit = polyval(coeffs, (xFit-UU(ii,i-1,1))./deltaX(ii));
        figure(1);hold on;
        plot(xFit, pFit, 'b-');
        drawnow;
    end
    
end

figure(1);
saveas(gcf, 'algebraic_mesh_higherRes', 'pdf');