clc;
clear
close all;

%% Control Parameters
tol = 1e-4;
omega = 1; % change later when norms are easier to find?
cBC = 0.5;

%% Initialize Computational Grid

% Computational Grid
nZeta = 51; % dZ = 1
nEta = 51; % dE = 1
[ZZ, EE] = meshgrid(linspace(0,1,nZeta), linspace(0,1,nEta));
dZ = 1/(nZeta-1);
dE = 1/(nEta-1);

%% Draw Lines
% Physical Grid -> boundary conditions
bAxis = 1;
aAxis = 2;
rEllipse = @(theta) aAxis*bAxis./sqrt((bAxis.*(cos(theta))).^2 + (aAxis.*sin(theta)).^2);
rMax = 5;
xRange = [0.5, 1];
xS = rEllipse(linspace(0,pi,nZeta)).*cos(linspace(0,pi,nZeta));
xN = rMax.*cos(linspace(0,pi,nZeta));
xW = linspace(xS(1), xN(1), nEta);
xE = linspace(xS(end), xN(end), nEta);

% yS = zeros(size(xS));
yN = rMax.*sin(linspace(0,pi,nZeta));
yS = rEllipse(linspace(0,pi,nZeta)).*sin(linspace(0,pi,nZeta));
% yN = ones(size(xN));
yW = linspace(yS(1), yN(1), nEta);
yE = linspace(yS(end), yN(end), nEta);

%% Build Physical Grid

UU = cat(3, ZZ, EE, zeros(size(ZZ)), zeros(size(EE))); % [X, Y, P, Q]
for i = 1:size(UU,2) % march west to east, intialize X and Y
   UU(:,i,1:2) = cat(3, linspace(xS(i), xN(i), nEta), linspace(yS(i), yN(i), nEta));
end

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

% figure(1);hold on;
% quiver(xN',yN',normU(:,1), normU(:,2)); hold on;
% quiver(xS',yS',normL(:,1), normL(:,2));
% 
% close all;
figure(1);
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal;

%% Fit P and Q

for i = 2:(size(UU,2)-1) % march from west to east, initialize P and Q
    % Get Derivatives
    xZ = 0.5.*(UU(2:end-1, i+1, 1) - UU(2:end-1, i-1,1))./dZ; % x_z
    yZ = 0.5.*(UU(2:end-1, i+1, 2) - UU(2:end-1, i-1,2))./dZ; % y_z
    xE = 0.5.*(UU(3:end,i,1) - UU(1:end-2,i,1))./dE; % x_e
    yE = 0.5.*(UU(3:end,i,2) - UU(1:end-2,i,2))./dE; % y_e

    % Calculate alpha, beta, gamma, jacobian
    alpha = xE.^2 + yE.^2;
    beta = xZ.*xE + yZ.*yE;
    gamma = xZ.^2 + yZ.^2;
    jac = xZ.*yE - yZ.*xE;
    
    % Calculate P and Q
    D_xy = -(alpha.*(UU(2:end-1, i+1,1:2)-2.*UU(2:end-1,i,1:2)+UU(2:end-1,i-1,1:2))./(dZ^2)...
            - 2.*beta.*0.5.*(0.5.*(UU(3:end,i+1,1:2) - UU(3:end,i-1,1:2))./dZ - 0.5.*(UU(1:end-2,i+1,1:2) - UU(1:end-2,i-1,1:2))./dZ)./dE ...
            + gamma.*(UU(3:end,i,1:2)-2.*UU(2:end-1,i,1:2)+UU(1:end-2,i,1:2))./(dE^2))./(jac.^2);
    UU(2:end-1,i,3:4) = (cat(3,yE,xZ).*D_xy - cat(3,xE,yZ).*D_xy(:,:,2:-1:1))./jac;
end

% % Pre-iterated Grid
figure(1);
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
figure(2);
plot(UU(:,:,3), UU(:,:,4), 'b-', UU(:,:,3)', UU(:,:,4)', 'b-'); axis equal; 
drawnow;

%% Line SOR
res = 1;
corr = 1;

while res(end) > tol && corr(end) > tol
    UUold = UU(:,:,1:2);
    %% Sweep through field 
%     if mod(length(res),2)
%         iSweep = 2:(size(UU,2)-1);
%     else
%         iSweep = fliplr(2:(size(UU,2)-1));
%     end
    iSweep = 2:(size(UU,2)-1);
    
    for i = iSweep
               
        % Update Boundary Conditions for p and q -> use interior points
        if (i == 2) || (i == (size(UU,2)-1))
            if i == 2
                ind = i-1;
%                 dZdx_dZdy = [(-1.5.*ZZ(2:end-1,1)+2.*ZZ(2:end-1,2)-0.5.*ZZ(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)),... dZ/dx
%                             (ZZ(3:end,1) - ZZ(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2))]; % dZ/dy
                        
                if any(abs(UU(3:end,1,2)-UU(1:end-2,1,2))< eps)
                    if any(abs(EE(3:end,1)-EE(1:end-2,1))< eps) 
                        dEdx_dEdy = [(-1.5.*EE(2:end-1,1)+2.*EE(2:end-1,2)-0.5.*EE(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)),... dE/dx
                                    zeros(size(EE(2:end-1,1)))]; % dE/dy
                    else
                       error('MESH NOT COMPATIBLE\n'); 
                    end
                elseif any(abs(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1))< eps) 
                    if any(abs(-1.5.*EE(2:end-1,1)+2.*EE(2:end-1,2)-0.5.*EE(2:end-1,3))< eps)
                        dEdx_dEdy = [zeros(size(EE(2:end-1,1))),... dE/dx
                                    (EE(3:end,1) - EE(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2))]; % dE/dy
                    else
                       error('MESH NOT COMPATIBLE\n'); 
                    end
                    
                else
                    dEdx_dEdy = [(-1.5.*EE(2:end-1,1)+2.*EE(2:end-1,2)-0.5.*EE(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)),... dE/dx
                                (EE(3:end,1) - EE(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2))]; % dE/dy
                end
                        
                normal = [(-(UU(3:end,1,2)-UU(1:end-2,1,2))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2)),... % n_x
                            (UU(3:end,1,1)-UU(1:end-2,1,1))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2)]; % n_y 
                        
                xZeta_WE =(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1))./dZ;
                yZeta_WE =(-1.5.*UU(2:end-1,1,2)+2.*UU(2:end-1,2,2)-0.5.*UU(2:end-1,3,2))./dZ;
                xEta_WE = 0.5.*(UU(3:end,ind,1) - UU(1:end-2,ind,1))./dE;
                yEta_WE = 0.5.*(UU(3:end,ind,2) - UU(1:end-2,ind,2))./dE;
                
                alpha_WE = xEta_WE.^2 + yEta_WE.^2;
                beta_WE = xZeta_WE.*xEta_WE + yZeta_WE.*yEta_WE;
                gamma_WE = xZeta_WE.^2 + yZeta_WE.^2;
                jac_WE = xZeta_WE.*yEta_WE - yZeta_WE.*xEta_WE;
                
                 % -(ax_zz - 2Bx_ze + gx_ee)/J^2
                D_xyZ = -( alpha_WE.*((UU(2:end-1,1,1:2)-2.*UU(2:end-1,2,1:2)+UU(2:end-1,3,1:2)))./(dZ^2)...
                    -2.*beta_WE.*((-1.5.*UU(3:end,1,1:2)+2.*UU(3:end,2,1:2)-0.5.*UU(3:end,3,1:2))./dZ-(-1.5.*UU(1:end-2,1,1:2)+2.*UU(1:end-2,2,1:2)-0.5.*UU(1:end-2,3,1:2))./dZ)./(2*dE) ...
                    +gamma_WE.*((UU(3:end,1,1:2)-2.*UU(2:end-1,1,1:2)+UU(1:end-2,1,1:2)))./(dE^2)  )./(jac_WE.^2);
            
            elseif i == (size(UU,2)-1)
                ind = i+1;             
%                 dZdx_dZdy = [(1.5.*ZZ(2:end-1,end)-2.*ZZ(2:end-1,end-1)+0.5.*ZZ(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)),... dZ/dx
%                     (ZZ(3:end,end) - ZZ(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2))]; % dZ/dy  

                if any(abs(UU(3:end,end,2)-UU(1:end-2,end,2))< eps)
                    if any((EE(3:end,end) - EE(1:end-2,end))< eps)
                        dEdx_dEdy = [(1.5.*EE(2:end-1,end)-2.*EE(2:end-1,end-1)+0.5.*EE(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)),... dE/dx
                                    zeros(size(EE(2:end-1,1)))]; % dE/dy
                    else
                       error('MESH NOT COMPATIBLE\n'); 
                    end
                elseif any(abs(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1))< eps)
                    if any(abs(1.5.*EE(2:end-1,end)-2.*EE(2:end-1,end-1)+0.5.*EE(2:end-1,end-2))< eps)
                        dEdx_dEdy = [zeros(size(EE(2:end-1,1))),... dE/dx
                                    (EE(3:end,end) - EE(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2))]; % dE/dy
                    else
                       error('MESH NOT COMPATIBLE\n'); 
                    end
                    
                else
                    dEdx_dEdy = [(1.5.*EE(2:end-1,end)-2.*EE(2:end-1,end-1)+0.5.*EE(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)),... dE/dx
                                    (EE(3:end,end) - EE(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2))]; % dE/dy
                end
                
                normal = [(-(UU(3:end,end,2)-UU(1:end-2,end,2))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2)),... % n_x
                            (UU(3:end,end,1)-UU(1:end-2,end,1))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2)]; % n_y
%                 normal = -[-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2)).^2)];
%                 normal = normal(2:end-1,:);
                        
                xZeta_WE =-(-1.5.*UU(2:end-1,end,1)+2.*UU(2:end-1,end-1,1)-0.5.*UU(2:end-1,end-2,1))./dZ;
                yZeta_WE =-(-1.5.*UU(2:end-1,end,2)+2.*UU(2:end-1,end-1,2)-0.5.*UU(2:end-1,end-2,2))./dZ;
                xEta_WE = 0.5.*(UU(3:end,ind,1) - UU(1:end-2,ind,1))./dE;
                yEta_WE = 0.5.*(UU(3:end,ind,2) - UU(1:end-2,ind,2))./dE;
                
                alpha_WE = xEta_WE.^2 + yEta_WE.^2;
                beta_WE = xZeta_WE.*xEta_WE + yZeta_WE.*yEta_WE;
                gamma_WE = xZeta_WE.^2 + yZeta_WE.^2;
                jac_WE = xZeta_WE.*yEta_WE - yZeta_WE.*xEta_WE;
                
                 % -(ax_zz - 2Bx_ze + gx_ee)/J^2
                D_xyZ = -( alpha_WE.*((UU(2:end-1,end,1:2)-2.*UU(2:end-1,end-1,1:2)+UU(2:end-1,end-2,1:2)))./(dZ^2)...
                    -2.*beta_WE.*(-(-1.5.*UU(3:end,end,1:2)+2.*UU(3:end,end-1,1:2)-0.5.*UU(3:end,end-2,1:2))./dZ+(-1.5.*UU(1:end-2,end,1:2)+2.*UU(1:end-2,end-1,1:2)-0.5.*UU(1:end-2,end-2,1:2))./dZ)./(2*dE) ...
                    +gamma_WE.*((UU(3:end,end,1:2)-2.*UU(2:end-1,end,1:2)+UU(1:end-2,end,1:2)))./(dE^2)  )./(jac_WE.^2);
            end

            % East/West Boundary
            UU(2:end-1, ind,3) = (yEta_WE.*D_xyZ(:,:,1) - xEta_WE.*D_xyZ(:,:,2))./jac_WE;% - cBC.*(sum(normal.*dZdx_dZdy,2));
            UU(2:end-1, ind,4) = (xZeta_WE.*D_xyZ(:,:,2) - yZeta_WE.*D_xyZ(:,:,1))./jac_WE - cBC.*(sum(normal.*dEdx_dEdy,2));
        end
        
        % Get Derivatives?
        xZ = 0.5.*(UU(2:end-1, i+1, 1) - UU(2:end-1, i-1,1))./dZ; % x_z
        yZ = 0.5.*(UU(2:end-1, i+1, 2) - UU(2:end-1, i-1,2))./dZ; % y_z
        xE = 0.5.*(UU(3:end,i,1) - UU(1:end-2,i,1))./dE; % x_e
        yE = 0.5.*(UU(3:end,i,2) - UU(1:end-2,i,2))./dE; % y_e
        UU_ZE = (diff(UU(3:end,[i+1,i-1],:), 1, 2)-diff(UU(1:end-2,[i+1,i-1],:), 1, 2))./(4*dZ*dE);

        % Calculate alpha, beta, gamma, jacobian
        alpha = xE.^2 + yE.^2;
        beta = xZ.*xE + yZ.*yE;
        gamma = xZ.^2 + yZ.^2;
        jac = xZ.*yE - yZ.*xE;
            
        D_xyN = -( alpha(end).*((UU(end,i+1,1:2)-2.*UU(end,i,1:2)+UU(end,i-1,1:2)))./(dZ^2)...
                -2.*beta(end).*UU_ZE(end,:,1:2) ...
                +gamma(end).*((UU(end,i,1:2)-2.*UU(end-1,i,1:2)+UU(end-2,i,1:2)))./(dE^2)  )./(jac(end)^2);

        D_xyS = -( alpha(1).*((UU(1,i+1,1:2)-2.*UU(1,i,1:2)+UU(1,i-1,1:2)))./(dZ^2)...
                -2.*beta(1).*UU_ZE(1,:,1:2) ...
                +gamma(1).*((UU(3,i,1:2)-2.*UU(2,i,1:2)+UU(1,i,1:2)))./(dE^2)  )./(jac(1)^2);

        % calculate Normals
        normalN = [(-(UU(end,i+1,2)-UU(end,i-1,2))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1))^2)),...
                    ((UU(end,i+1,1)-UU(end,i-1,1))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1)).^2))];
        normalS = [(-(UU(1,i+1,2)-UU(1,i-1,2))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1))^2)),...
                    ((UU(1,i+1,1)-UU(1,i-1,1))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1)).^2))];

        % North Boundary
        dZdx_dZdy_N = [(ZZ(end,i+1) - ZZ(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1)), (1.5*ZZ(end,i)-2*ZZ(end-1,i)+0.5*ZZ(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2))];
        dEdx_dEdy_N = [(EE(end,i+1) - EE(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1)), (1.5*EE(end,i)-2*EE(end-1,i)+0.5*EE(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2))];
        UU(end,i,3) = (yE(end)*D_xyN(:,:,1) - xE(end)*D_xyN(:,:,2))./jac(end) - cBC.*(sum(normalN.*dZdx_dZdy_N));
        UU(end,i,4) = (xZ(end)*D_xyN(:,:,2) - yZ(end)*D_xyN(:,:,1))./jac(end);% - cBC.*(sum(normalN.*dEdx_dEdy_N));

        % South Boundary
        dZdx_dZdy_S = [(ZZ(1,i+1) - ZZ(1,i-1))./(UU(1,i+1,1) - UU(1,i-1,1)), (1.5*ZZ(1,i)-2*ZZ(2,i)+0.5*ZZ(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2))];
        dEdx_dEdy_S = [(EE(1,i+1) - EE(1,i-1))./(UU(1,i+1,1) - UU(1,i-1,1)),  (1.5*EE(1,i)-2*EE(2,i)+0.5*EE(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2))];
        UU(1,i,3) = (yE(1)*D_xyS(:,:,1) - xE(1)*D_xyS(:,:,2))./jac(1) - cBC.*(sum(normalS.*dZdx_dZdy_S));
        UU(1,i,4) = (xZ(1)*D_xyS(:,:,2) - yZ(1)*D_xyS(:,:,1))./jac(1);% - cBC.*(sum(normalS.*dEdx_dEdy_S));

        % Calculate thomas coefficients
        aa = repmat([0; gamma; 0],4,1)./(dE^2);
        bb = repmat([1; -2.*(alpha./(dZ^2) + gamma./(dE^2)); 1], size(UU,3),1);
        cc = repmat([0; gamma; 0],4,1)./(dE^2);
        dd = [  UU(1,i,1); -alpha.*sum(UU(2:end-1,[i+1,i-1],1), 2)./(dZ^2) + 2.*beta.*UU_ZE(:,:,1) - (jac.^2).*(UU(2:end-1,i,3).*xZ + UU(2:end-1,i,4).*xE); UU(end,i,1);...
                UU(1,i,2); -alpha.*sum(UU(2:end-1,[i+1,i-1],2), 2)./(dZ^2) + 2.*beta.*UU_ZE(:,:,2) - (jac.^2).*(UU(2:end-1,i,3).*yZ + UU(2:end-1,i,4).*yE); UU(end,i,2);...
                UU(1,i,3); -alpha.*sum(UU(2:end-1,[i+1,i-1],3), 2)./(dZ^2) + 2.*beta.*UU_ZE(:,:,3); UU(end,i,3);...
                UU(1,i,4); -alpha.*sum(UU(2:end-1,[i+1,i-1],4), 2)./(dZ^2) + 2.*beta.*UU_ZE(:,:,4); UU(end,i,4)];
            
        UU(:,i,:) = omega.*reshape(thomas3(aa, bb, cc, dd),size(UU,1),1,size(UU,3)) + (1-omega).*UU(:,i,:);
        
    end
    
    %% Calculate Residual
    res(end+1) = meshResidual(UU);
    corr(end+1) = max(max(max(abs(UU(:,:,1:2) - UUold),[],3)));
    
    figure(1);
    plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
    hold on;
    plot(UU([1 size(UU,1)],:,1)', UU([1 size(UU,1)],:,2)', 'b-'); axis equal;
%     normExit = -[-1./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2)).^2), polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat(1:end-1,end),UU(:,end,2)).^2)];
%     quiver(xN',yN', normU(:,1), normU(:,2));
%     quiver(UU(:,end,1), UU(:,end,2), normExit(:,1), normExit(:,2)); hold off;
    
%     figure(2);
%     plot(UU(:,:,3), UU(:,:,4), 'b-', UU(:,:,3)', UU(:,:,4)', 'b-'); axis equal; 
    hold off;
    drawnow;
end

figure();
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
title('Biharmonic Mesh -- Nozzle');
saveas(gcf, 'biharmonic_mesh', 'pdf');