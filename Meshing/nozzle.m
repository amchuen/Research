clc;
clear
close all;

%% Control Parameters
tol = 1e-5;
omega = 1; % change later when norms are easier to find?
cBC = 1;

%% Initialize Computational Grid and Physical Grid

% Computational Grid
nZeta = 21; % dZ = 1
nEta = 7; % dE = 1
[ZZ, EE] = meshgrid(1:nZeta, 1:nEta);

% Physical Grid -> boundary conditions
xS = linspace(0.5, 2, nZeta);
xW = zeros(nEta, 1);
xN = xS;
xE = ones(nEta, 1);

yS = zeros(size(xS));
yN = [linspace(0.5,1, sum(xN<1)), ones(size(xN(xN >=1)))];%[0.5.*(1 + (2.*xN(xN <1)-1).^2), ones(size(xN(xN >=1)))];
yW = linspace(yS(1), yN(1), nEta);
yE = linspace(yS(end), yN(end), nEta);

UU = cat(3, ZZ, EE, zeros(size(ZZ)), zeros(size(EE))); % [X, Y, P, Q]
for i = 1:size(UU,2)
   UU(:,i,1:2) = cat(3, linspace(xS(i), xN(i), nEta), linspace(yS(i), yN(i), nEta));
end
for i = 2:(size(UU,2)-1) % march from west to east
    % Get Derivatives?
    xZeta = UU(2:end-1, i+1, 1) - UU(2:end-1, i-1,1);
    yZeta = UU(2:end-1, i+1, 2) - UU(2:end-1, i-1,2);
    xEta = UU(3:end,i,1) - UU(1:end-2,i,1);
    yEta = UU(3:end,i,2) - UU(1:end-2,i,2);

    % Calculate alpha, beta, gamma, jacobian
    alpha = xEta.^2 + yEta.^2;
    beta = xZeta.*xEta + yZeta.*yEta;
    gamma = xZeta.^2 + yZeta.^2;
    jac = xZeta.*yEta - yZeta.*xEta;
    
    if (i == 2) || (i == (size(UU,2)-1))
        if i == 2
            ind = i-1;
            dZdx = (-1.5.*ZZ(2:end-1,1)+2.*ZZ(2:end-1,2)-0.5.*ZZ(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)).*(-(UU(3:end,1,2)-UU(1:end-2,1,2))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2))...
                    + (ZZ(3:end,1) - ZZ(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2)).*(UU(3:end,1,1)-UU(1:end-2,1,1))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2); % normal is <-1,0>
            dEdx = (-1.5.*EE(2:end-1,1)+2.*EE(2:end-1,2)-0.5.*EE(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)).*(-(UU(3:end,1,2)-UU(1:end-2,1,2))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2))...
                    + (EE(3:end,1) - EE(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2)).*(UU(3:end,1,1)-UU(1:end-2,1,1))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2); % normal is <-1,0>
            alpha_WE = diff(UU(2:end-1,1:2,1),1,2).^2 + diff(UU(2:end-1,1:2,2),1,2).^2;
        elseif i == (size(UU,2)-1)
            ind = i+1;                
            dZdx = (1.5.*ZZ(2:end-1,end)-2.*ZZ(2:end-1,end-1)+0.5.*ZZ(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)).*(-(UU(3:end,end,2)-UU(1:end-2,end,2))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2))...
                    + (ZZ(3:end,end) - ZZ(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2)).*(UU(3:end,end,1)-UU(1:end-2,end,1))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2); % normal is <-1,0>
            dEdx = (1.5.*EE(2:end-1,end)-2.*EE(2:end-1,end-1)+0.5.*EE(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)).*(-(UU(3:end,end,2)-UU(1:end-2,end,2))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2))...
                    + (EE(3:end,end) - EE(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2)).*(UU(3:end,end,1)-UU(1:end-2,end,1))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2); % normal is <-1,0>
            alpha_WE = diff(UU(2:end-1,end-1:end,1),1,2).^2 + diff(UU(2:end-1,end-1:end,2),1,2).^2;
        end

        % -(ax_zz - 2Bx_ze + gx_ee)/J^2
        DxZ = -( alpha.*((UU(2:end-1,i+1,1)-2.*UU(2:end-1,i,1)+UU(2:end-1,i-1,1)))...
                -2.*beta.*(diff(UU(3:end,[i+1,i-1],1), 1, 2)-diff(UU(1:end-2,[i+1,i-1],1), 1, 2))...
                +gamma.*((UU(3:end,i,1)-2.*UU(2:end-1,i,1)+UU(1:end-2,i,1)))  )./(jac.^2);

        DyZ = -( alpha.*((UU(2:end-1,i+1,2)-2.*UU(2:end-1,i,2)+UU(2:end-1,i-1,2)))...
                -2.*beta.*(diff(UU(3:end,[i+1,i-1],2), 1, 2)-diff(UU(1:end-2,[i+1,i-1],2), 1, 2))...
                +gamma.*((UU(3:end,i,2)-2.*UU(2:end-1,i,2)+UU(1:end-2,i,2)))  )./(jac.^2);

        % East/West Boundary
        UU(2:end-1, ind,3) = (yEta.*DxZ - xEta.*DyZ)./jac;% - cBC.*(dZdx);
        UU(2:end-1, ind,4) = (xZeta.*DyZ - yZeta.*DxZ)./jac - cBC.*(dEdx);
        
        UU(:,ind,3:4) = cat(3, linspace(UU(1,ind,3),UU(end,ind,3),nEta), linspace(UU(1,ind,4),UU(end,ind,4),nEta));    
    end
    
    DxN = -( alpha(end).*((UU(end-1,i+1,1)-2.*UU(end-1,i,1)+UU(end-1,i-1,1)))...
            -2.*beta(end).*(diff(UU(end,[i+1,i-1],1), 1, 2)-diff(UU(end-2,[i+1,i-1],1), 1, 2))...
            +gamma(end).*((UU(end,i,1)-2.*UU(end-1,i,1)+UU(end-2,i,1)))  )./(jac(end)^2);

    DyN = -( alpha(end).*((UU(end-1,i+1,2)-2.*UU(end-1,i,2)+UU(end-1,i-1,2)))...
            -2.*beta(end).*(diff(UU(end,[i+1,i-1],2), 1, 2)-diff(UU(end-2,[i+1,i-1],2), 1, 2))...
            +gamma(end).*((UU(end,i,2)-2.*UU(end-1,i,2)+UU(end-2,i,2)))  )./(jac(end)^2);

    DxS = -( alpha(1).*((UU(2,i+1,1)-2.*UU(2,i,1)+UU(2,i-1,1)))...
            -2.*beta(1).*(diff(UU(3,[i+1,i-1],1), 1, 2)-diff(UU(1,[i+1,i-1],1), 1, 2))...
            +gamma(1).*((UU(3,i,1)-2.*UU(2,i,1)+UU(1,i,1)))  )./(jac(1)^2);

    DyS = -( alpha(1).*((UU(2,i+1,2)-2.*UU(2,i,2)+UU(2,i-1,2)))...
            -2.*beta(1).*(diff(UU(3,[i+1,i-1],2), 1, 2)-diff(UU(1,[i+1,i-1],2), 1, 2))...
            +gamma(1).*((UU(3,i,2)-2.*UU(2,i,2)+UU(1,i,2)))  )./(jac(1)^2);

    UU(end,i,3) = (yEta(end)*DxN - xEta(end)*DyN)./jac(end)...
                    - cBC.*((ZZ(end,i+1) - ZZ(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-(UU(end,i+1,2)-UU(end,i-1,2))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1))^2)) + (1.5*ZZ(end,i)-2*ZZ(end-1,i)+0.5*ZZ(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * ((UU(end,i+1,1)-UU(end,i-1,1))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1)).^2)));
    UU(end,i,4) = (xZeta(end)*DyN - yZeta(end)*DxN)./jac(end);%...
                    %- cBC.*((EE(end,i+1) - EE(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-(UU(end,i+1,2)-UU(end,i-1,2))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1))^2)) + (1.5*EE(end,i)-2*EE(end-1,i)+0.5*EE(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * ((UU(end,i+1,1)-UU(end,i-1,1))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1)).^2)));

%     UU(1,i,3) = (yEta(1)*DxS - xEta(1)*DyS)./jac(1) - cBC.*((-1.5*ZZ(1,i)+2*ZZ(2,i)-0.5*ZZ(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
%     UU(1,i,4) = (xZeta(1)*DyS - yZeta(1)*DxS)./jac(1) - cBC.*((-1.5*EE(1,i)+2*EE(2,i)-0.5*EE(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
    UU(1,i,3) = (yEta(1)*DxS - xEta(1)*DyS)./jac(1)...
                    - cBC.*((ZZ(1,i+1) - ZZ(1,i-1))./(UU(1,i+1,1) - UU(1,i-1,1))*(-(UU(1,i+1,2)-UU(1,i-1,2))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1))^2)) + (1.5*ZZ(1,i)-2*ZZ(2,i)+0.5*ZZ(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2)) * ((UU(1,i+1,1)-UU(1,i-1,1))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1)).^2)));
    UU(1,i,4) = (xZeta(1)*DyS - yZeta(1)*DxS)./jac(1);%...
                    %- cBC.*((EE(1,i+1) - EE(1,i-1))/(UU(1,i+1,1) - UU(1,i-1,1))*(-(UU(1,i+1,2)-UU(1,i-1,2))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1))^2)) + (1.5*EE(1,i)-2*EE(2,i)+0.5*EE(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2)) * ((UU(1,i+1,1)-UU(1,i-1,1))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1)).^2)));
    
    UU(:,i,3:4) = cat(3, linspace(UU(1,i,3),UU(end,i,3),nEta), linspace(UU(1,i,4),UU(end,i,4),nEta));    
end

% % Pre-iterated Grid
figure(1);
plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
figure(2);
plot(UU(:,:,3), UU(:,:,4), 'b-', UU(:,:,3)', UU(:,:,4)', 'b-'); axis equal; 
drawnow;

%% Line SOR
res = 1;

while res(end) > tol
    %% Sweep through field 
    if mod(length(res),2)
        iSweep = 2:(size(UU,2)-1);
    else
        iSweep = fliplr(2:(size(UU,2)-1));
    end
    
    for i = iSweep
        % Get Derivatives?
        xZeta = UU(2:end-1, i+1, 1) - UU(2:end-1, i-1,1);
        yZeta = UU(2:end-1, i+1, 2) - UU(2:end-1, i-1,2);
        xEta = UU(3:end,i,1) - UU(1:end-2,i,1);
        yEta = UU(3:end,i,2) - UU(1:end-2,i,2);
        
        % Calculate alpha, beta, gamma, jacobian
        alpha = xEta.^2 + yEta.^2;
        beta = xZeta.*xEta + yZeta.*yEta;
        gamma = xZeta.^2 + yZeta.^2;
        jac = xZeta.*yEta - yZeta.*xEta;
        
        % Update Boundary Conditions for p and q -> use interior points
        if (i == 2) || (i == (size(UU,2)-1))
            if i == 2
                ind = i-1;
                dZdx = (-1.5.*ZZ(2:end-1,1)+2.*ZZ(2:end-1,2)-0.5.*ZZ(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)).*(-(UU(3:end,1,2)-UU(1:end-2,1,2))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2))...
                        + (ZZ(3:end,1) - ZZ(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2)).*(UU(3:end,1,1)-UU(1:end-2,1,1))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2); % normal is <-1,0>
                dEdx = (-1.5.*EE(2:end-1,1)+2.*EE(2:end-1,2)-0.5.*EE(2:end-1,3))./(-1.5.*UU(2:end-1,1,1)+2.*UU(2:end-1,2,1)-0.5.*UU(2:end-1,3,1)).*(-(UU(3:end,1,2)-UU(1:end-2,1,2))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2))...
                        + (EE(3:end,1) - EE(1:end-2,1))./(UU(3:end,1,2)-UU(1:end-2,1,2)).*(UU(3:end,1,1)-UU(1:end-2,1,1))./sqrt((UU(3:end,1,2)-UU(1:end-2,1,2)).^2+(UU(3:end,1,1)-UU(1:end-2,1,1)).^2); % normal is <-1,0>
                alpha_WE = diff(UU(2:end-1,1:2,1),1,2).^2 + diff(UU(2:end-1,1:2,2),1,2).^2;
            elseif i == (size(UU,2)-1)
                ind = i+1;                
                dZdx = (1.5.*ZZ(2:end-1,end)-2.*ZZ(2:end-1,end-1)+0.5.*ZZ(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)).*(-(UU(3:end,end,2)-UU(1:end-2,end,2))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2))...
                        + (ZZ(3:end,end) - ZZ(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2)).*(UU(3:end,end,1)-UU(1:end-2,end,1))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2); % normal is <-1,0>
                dEdx = (1.5.*EE(2:end-1,end)-2.*EE(2:end-1,end-1)+0.5.*EE(2:end-1,end-2))./(1.5.*UU(2:end-1,end,1)-2.*UU(2:end-1,end-1,1)+0.5.*UU(2:end-1,end-2,1)).*(-(UU(3:end,end,2)-UU(1:end-2,end,2))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2))...
                        + (EE(3:end,end) - EE(1:end-2,end))./(UU(3:end,end,2)-UU(1:end-2,end,2)).*(UU(3:end,end,1)-UU(1:end-2,end,1))./sqrt((UU(3:end,end,2)-UU(1:end-2,end,2)).^2+(UU(3:end,end,1)-UU(1:end-2,end,1)).^2); % normal is <-1,0>
                alpha_WE = diff(UU(2:end-1,end-1:end,1),1,2).^2 + diff(UU(2:end-1,end-1:end,2),1,2).^2;
            end
            
            % -(ax_zz - 2Bx_ze + gx_ee)/J^2
            DxZ = -( alpha.*((UU(2:end-1,i+1,1)-2.*UU(2:end-1,i,1)+UU(2:end-1,i-1,1)))...
                    -2.*beta.*(diff(UU(3:end,[i+1,i-1],1), 1, 2)-diff(UU(1:end-2,[i+1,i-1],1), 1, 2))...
                    +gamma.*((UU(3:end,i,1)-2.*UU(2:end-1,i,1)+UU(1:end-2,i,1)))  )./(jac.^2);
                
            DyZ = -( alpha.*((UU(2:end-1,i+1,2)-2.*UU(2:end-1,i,2)+UU(2:end-1,i-1,2)))...
                    -2.*beta.*(diff(UU(3:end,[i+1,i-1],2), 1, 2)-diff(UU(1:end-2,[i+1,i-1],2), 1, 2))...
                    +gamma.*((UU(3:end,i,2)-2.*UU(2:end-1,i,2)+UU(1:end-2,i,2)))  )./(jac.^2);
            
            % East/West Boundary
            UU(2:end-1, ind,3) = (yEta.*DxZ - xEta.*DyZ)./jac;% - cBC.*(dZdx);
            UU(2:end-1, ind,4) = (xZeta.*DyZ - yZeta.*DxZ)./jac - cBC.*(dEdx);
        end
            
%             % North Boundary
%             UU(end,i,3) = (yEta(end)*DxZ(end) - xEta(end)*DyZ(end))./jac(end)...
%                             - cBC.*((ZZ(end,i+1) - ZZ(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-2*UU(end,i,1)./sqrt(1+4*UU(end,i,1)^2)) + (1.5*ZZ(end,i)-2*ZZ(end-1,i)+0.5*ZZ(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * (1./sqrt(1+4*UU(end,i,1)^2)));
%             UU(end,i,4) = (xZeta(end)*DyZ(end) - yZeta(end)*DxZ(end))./jac(end)...
%                             - cBC.*((EE(end,i+1) - EE(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-2*UU(end,i,1)./sqrt(1+4*UU(end,i,1)^2)) + (1.5*EE(end,i)-2*EE(end-1,i)+0.5*EE(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * (1./sqrt(1+4*UU(end,i,1)^2)));
%             
%             % South Boundary
%             UU(1,i,3) = (yEta(1)*DxZ(1) - xEta(1)*DyZ(1))./jac(1) - cBC.*((-1.5*ZZ(1,i)+2*ZZ(2,i)-0.5*ZZ(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
%             UU(1,i,4) = (xZeta(1)*DyZ(1) - yZeta(1)*DxZ(1))./jac(1) - cBC.*((-1.5*EE(1,i)+2*EE(2,i)-0.5*EE(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
        DxN = -( alpha(end).*((UU(end-1,i+1,1)-2.*UU(end-1,i,1)+UU(end-1,i-1,1)))...
                -2.*beta(end).*(diff(UU(end,[i+1,i-1],1), 1, 2)-diff(UU(end-2,[i+1,i-1],1), 1, 2))...
                +gamma(end).*((UU(end,i,1)-2.*UU(end-1,i,1)+UU(end-2,i,1)))  )./(jac(end)^2);

        DyN = -( alpha(end).*((UU(end-1,i+1,2)-2.*UU(end-1,i,2)+UU(end-1,i-1,2)))...
                -2.*beta(end).*(diff(UU(end,[i+1,i-1],2), 1, 2)-diff(UU(end-2,[i+1,i-1],2), 1, 2))...
                +gamma(end).*((UU(end,i,2)-2.*UU(end-1,i,2)+UU(end-2,i,2)))  )./(jac(end)^2);

        DxS = -( alpha(1).*((UU(2,i+1,1)-2.*UU(2,i,1)+UU(2,i-1,1)))...
                -2.*beta(1).*(diff(UU(3,[i+1,i-1],1), 1, 2)-diff(UU(1,[i+1,i-1],1), 1, 2))...
                +gamma(1).*((UU(3,i,1)-2.*UU(2,i,1)+UU(1,i,1)))  )./(jac(1)^2);

        DyS = -( alpha(1).*((UU(2,i+1,2)-2.*UU(2,i,2)+UU(2,i-1,2)))...
                -2.*beta(1).*(diff(UU(3,[i+1,i-1],2), 1, 2)-diff(UU(1,[i+1,i-1],2), 1, 2))...
                +gamma(1).*((UU(3,i,2)-2.*UU(2,i,2)+UU(1,i,2)))  )./(jac(1)^2);

        UU(end,i,3) = (yEta(end)*DxN - xEta(end)*DyN)./jac(end)...
                        - cBC.*((ZZ(end,i+1) - ZZ(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-(UU(end,i+1,2)-UU(end,i-1,2))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1))^2)) + (1.5*ZZ(end,i)-2*ZZ(end-1,i)+0.5*ZZ(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * ((UU(end,i+1,1)-UU(end,i-1,1))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1)).^2)));
        UU(end,i,4) = (xZeta(end)*DyN - yZeta(end)*DxN)./jac(end);...
                        %- cBC.*((EE(end,i+1) - EE(end,i-1))/(UU(end,i+1,1) - UU(end,i-1,1))*(-(UU(end,i+1,2)-UU(end,i-1,2))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1))^2)) + (1.5*EE(end,i)-2*EE(end-1,i)+0.5*EE(end-2,i))/(1.5*UU(end,i,2)-2*UU(end-1,i,2)+0.5*UU(end-2,i,2)) * ((UU(end,i+1,1)-UU(end,i-1,1))./sqrt((UU(end,i+1,2)-UU(end,i-1,2)).^2+(UU(end,i+1,1)-UU(end,i-1,1)).^2)));

%         UU(1,i,3) = (yEta(1)*DxS - xEta(1)*DyS)./jac(1) - cBC.*((-1.5*ZZ(1,i)+2*ZZ(2,i)-0.5*ZZ(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
%         UU(1,i,4) = (xZeta(1)*DyS - yZeta(1)*DxS)./jac(1) - cBC.*((-1.5*EE(1,i)+2*EE(2,i)-0.5*EE(3,i))./(-1.5*UU(1,i,2)+2*UU(2,i,2)-0.5*UU(3,i,2)));
        UU(1,i,3) = (yEta(1)*DxS - xEta(1)*DyS)./jac(1)...
                        - cBC.*((ZZ(1,i+1) - ZZ(1,i-1))./(UU(1,i+1,1) - UU(1,i-1,1))*(-(UU(1,i+1,2)-UU(1,i-1,2))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1))^2)) + (1.5*ZZ(1,i)-2*ZZ(2,i)+0.5*ZZ(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2)) * ((UU(1,i+1,1)-UU(1,i-1,1))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1)).^2)));
        UU(1,i,4) = (xZeta(1)*DyS - yZeta(1)*DxS)./jac(1);...
                        %- cBC.*((EE(1,i+1) - EE(1,i-1))/(UU(1,i+1,1) - UU(1,i-1,1))*(-(UU(1,i+1,2)-UU(1,i-1,2))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1))^2)) + (1.5*EE(1,i)-2*EE(2,i)+0.5*EE(3,i))/(1.5*UU(1,i,2)-2*UU(2,i,2)+0.5*UU(3,i,2)) * ((UU(1,i+1,1)-UU(1,i-1,1))./sqrt((UU(1,i+1,2)-UU(1,i-1,2)).^2+(UU(1,i+1,1)-UU(1,i-1,1)).^2)));

        % Calculate thomas coefficients
        aa = repmat([0; gamma; 0],4,1);
        bb = repmat([1; -2.*(alpha + gamma); 1], size(UU,3),1);
        cc = repmat([0; gamma; 0],4,1);
        dd = [  UU(1,i,1); -alpha.*sum(UU(2:end-1,[i+1,i-1],1), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],1), 1, 2)-diff(UU(1:end-2,[i+1,i-1],1), 1, 2)) - (jac.^2).*(UU(2:end-1,i,3).*xZeta + UU(2:end-1,i,4).*xEta); UU(end,i,1);...
                UU(1,i,2); -alpha.*sum(UU(2:end-1,[i+1,i-1],2), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],2), 1, 2)-diff(UU(1:end-2,[i+1,i-1],2), 1, 2)) - (jac.^2).*(UU(2:end-1,i,3).*yZeta + UU(2:end-1,i,4).*yEta); UU(end,i,2);...
                UU(1,i,3); -alpha.*sum(UU(2:end-1,[i+1,i-1],3), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],3), 1, 2)-diff(UU(1:end-2,[i+1,i-1],3), 1, 2)); UU(end,i,3);...
                UU(1,i,4); -alpha.*sum(UU(2:end-1,[i+1,i-1],4), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],4), 1, 2)-diff(UU(1:end-2,[i+1,i-1],4), 1, 2)); UU(end,i,4)];
            
        UU(:,i,:) = omega.*reshape(thomas3(aa, bb, cc, dd),size(UU,1),1,size(UU,3)) + (1-omega).*UU(:,i,:);
        
    end
    
    %% Calculate Residual
    res(end+1) = meshResidual(UU);
    
    figure(1);
    plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal; 
    figure(2);
    plot(UU(:,:,3), UU(:,:,4), 'b-', UU(:,:,3)', UU(:,:,4)', 'b-'); axis equal; 
    drawnow;
end