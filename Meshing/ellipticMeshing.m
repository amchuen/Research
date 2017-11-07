clc;
clear;
close all;

%% Generate Computational Grid

nE = 21; % number of Eta points -> y-dirs
nZ = 101; % number of Zeta points -> x dir

de = 1./(nE-1);
dz = 1./(nZ-1);

[ZZ, EE] = meshgrid(0:dz:1, 0:de:1);

%% Generate Domain and Boundary Conditions

XX = ones(nE, nZ);
YY = ones(nE, nZ);
rad = 1;
tau = rad;
betaM = 1.05; % mesh clustering

% E1 - west
XX(:,1) = -5;
YY(:,1) = linspace(0,2.5, nE);

% E2 - east
XX(:,end) = 5;
YY(:,end) = YY(:,1);

% E3 - south
XX(1,:) = linspace(XX(1,1), XX(1,end), nZ);
YY(1,:) = [ zeros(size(XX(1,XX(1,:) < -rad))),...
            tau.*sqrt(rad^2 - XX(1,XX(1,:) >= -rad & XX(1,:) <= rad).^2),...
            zeros(size(XX(1,XX(1,:) > rad)))];

% E4 - north
XX(end,:) = XX(1,:);
YY(end,:) = 2.5;

% Fill in rest of values
for i = 1:(size(XX,2))
   XX(:,i) = linspace(XX(1,i),XX(end,i),length(XX(:,i))); 
   HH = YY(end,i) - YY(1,i);
   YY(2:end-1,i) = YY(1,i) + HH .* ((betaM+1)-(betaM-1).*(((betaM+1)/(betaM-1)).^(1-EE(2:end-1,i))))./(((betaM+1)./(betaM-1)).^(1-EE(2:end-1,i))+1);
%    YY(:,i) = linspace(YY(1,i),YY(end,i),length(YY(:,i)));
    
end

figure();plot(XX, YY, 'b-', XX', YY', 'b-');axis equal
res = 1;
tol = 1e-3;
wSOR = 1.74;

% XXold = XX;
% YYold = YY;

%% Compute Mesh with Line SOR?

while res(end) > tol
    XXold = XX;
    YYold = YY;
% XX = XXold;
%     YY = YYold;
    
    for i = 2:(size(XX,2)-1)
        xE = [0; (XX(3:end,i) - XX(1:end-2,i) )./(2*de); 0];
        xZ = [0; (XX(2:end-1,i+1) - XX(2:end-1,i-1))./(2*dz); 0];
        yE = [0; (YY(3:end,i) - YY(1:end-2,i) )./(2*de); 0];
        yZ = [0; (YY(3:end,i+1) - YY(1:end-2,i-1))./(2*de); 0];
        xZE = [0; (XX(3:end,i+1)-XX(1:end-2,i+1)-XX(3:end,i-1)+XX(1:end-2,i-1))./(4*de*dz);0];
        yZE = [0; (YY(3:end,i+1)-YY(1:end-2,i+1)-YY(3:end,i-1)+YY(1:end-2,i-1))./(4*de*dz);0];

        eXX = [0;((EE(2:end-1,i+1) - EE(2:end-1,i))./(XX(2:end-1,i+1) - XX(2:end-1,i)) - (EE(2:end-1,i) - EE(2:end-1,i-1))./(XX(2:end-1,i) - XX(2:end-1,i-1)))./(0.5.*(XX(2:end-1,i+1) - XX(2:end-1,i-1))); 0];
        eYY = [0; ((EE(3:end,i) - EE(2:end-1,i))./(YY(3:end,i) - YY(2:end-1,i)) - (EE(2:end-1,i) - EE(1:end-2,i))./(YY(2:end-1,i) - YY(1:end-2,i)))./(0.5.*(YY(3:end,i) - YY(1:end-2,i))); 0];
        zXX = [0; ((ZZ(2:end-1,i+1) - ZZ(2:end-1,i))./(XX(2:end-1,i+1) - XX(2:end-1,i)) - (ZZ(2:end-1,i) - ZZ(2:end-1,i-1))./(XX(2:end-1,i) - XX(2:end-1,i-1)))./(0.5.*(XX(2:end-1,i+1) - XX(2:end-1,i-1))); 0];
        zYY = [0; ((ZZ(3:end,i) - ZZ(2:end-1,i))./(YY(3:end,i) - YY(2:end-1,i)) - (ZZ(2:end-1,i) - ZZ(1:end-2,i))./(YY(2:end-1,i) - YY(1:end-2,i)))./(0.5.*(YY(3:end,i) - YY(1:end-2,i))); 0];

        alpha = xE.^2 + yE.^2;
        beta = xZ.*xE + yZ.*yE;
        gamma = xZ.^2 + yZ.^2;
        Jac = xZ.*yE - xE.*yZ;

        % Calculate the Laplacian using Tridiagonal Matrix
        aa = gamma ./ (de^2);
        bb = [1; -2.*(alpha(2:end-1)./(dz^2) + gamma(2:end-1)./(de^2)); 1];
        cc = aa;
%         - Jac(2:end-1).^2.*((eXX(2:end-1)+eYY(2:end-1)).*xZ(2:end-1) + (zXX(2:end-1)+zYY(2:end-1)).*xE(2:end-1))
        ddX = [XX(1,i); 2.*beta(2:end-1).*xZE(2:end-1) - alpha(2:end-1)./(dz^2).*(XX(2:end-1,i-1)+XX(2:end-1,i+1)); XX(end,i)];
        ddY = [YY(1,i); 2.*beta(2:end-1).*yZE(2:end-1) - alpha(2:end-1)./(dz^2).*(YY(2:end-1,i-1)+YY(2:end-1,i+1)); YY(end,i)];

        XX(:,i) = wSOR.*thomas3(aa,bb,cc,ddX) + (1-wSOR).*XX(:,i);
        YY(:,i) = wSOR.*thomas3(aa,bb,cc,ddY) + (1-wSOR).*YY(:,i);
        
        resCheck = max(max(abs(XX(:,i)-XXold(:,i))), max(abs(YY(:,i)-YYold(:,i))));
        
        if resCheck > 1
           error('Residual is growing too fast?\n'); 
        end

    end
    
    res(end+1) = max(max(abs(XX(:)-XXold(:))), max(abs(YY(:)-YYold(:))));
%     if res(end) > 1
%        error('Residual is growing too fast?\n'); 
%         
%     end
end


%% Post-Process?

figure();plot(XX, YY, 'b-', XX', YY', 'b-');axis equal

figure();semilogy(1:length(res), res);