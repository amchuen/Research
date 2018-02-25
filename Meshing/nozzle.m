clc;
clear
close all;

%% Initialize Computational Grid and Physical Grid

% Computational Grid
nZeta = 101; % dZ = 1
nEta = 51; % dE = 1
[ZZ, EE] = meshgrid(1:nZeta, 1:nEta);

% Physical Grid -> boundary conditions
xS = linspace(0, 1, nZeta);
xW = zeros(nEta, 1);
xN = xS;
xE = ones(nEta, 1);

yS = zeros(size(xS));
yN = 0.5.*(1 + (2.*xN-1).^2);
yW = linspace(yS(1), yN(1), nEta);
yE = linspace(yS(end), yN(end), nEta);

UU = cat(3, ZZ, EE, zeros(size(ZZ)), zeros(size(EE))); % [X, Y, P, Q]
for i = 1:size(UU,2) % march from west to east
    UU(:,i,1:2) = cat(3, linspace(xS(i), xN(i), nEta), linspace(yS(i), yN(i), nEta));    
end

% % Pre-iterated Grid
% figure();
% plot(UU(:,:,1), UU(:,:,2), 'b-', UU(:,:,1)', UU(:,:,2)', 'b-'); axis equal;

%% Control Parameters
tol = 1e-5;
omega = 1; % change later when norms are easier to find?
cBC = 1;

%% Line SOR
res = 1;

while res(end) > tol
    %% Sweep through field 
    for i = 2:(size(UU,2)-1)
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
            DxZ = -( alpha.*((UU(2:end-1,i+1,1)-2.*UU(2:end-1,i,1)+UU(2:end-1,i-1,1)))...
                    -2.*beta.*(diff(UU(3:end,[i+1,i-1],1), 1, 2)-diff(UU(1:end-2,[i+1,i-1],1), 1, 2))...
                    +gamma.*((UU(3:end,i,1)-2.*UU(2:end-1,i,1)+UU(1:end-2,i,1)))  )./(jac.^2);
                
            DyZ = -( alpha.*((UU(2:end-1,i+1,2)-2.*UU(2:end-1,i,2)+UU(2:end-1,i-1,2)))...
                    -2.*beta.*(diff(UU(3:end,[i+1,i-1],2), 1, 2)-diff(UU(1:end-2,[i+1,i-1],2), 1, 2))...
                    +gamma.*((UU(3:end,i,2)-2.*UU(2:end-1,i,2)+UU(1:end-2,i,2)))  )./(jac.^2);
            
            if i == 2
                ind = i-1;
%                 dZdx = 1./diff(UU(2:end-1,1:2,1),1,2);
%                 dEdy = 1./(UU(3:end,1,2) - UU(1:end-2,1,2));
            elseif i == (size(UU,2)-1)
                ind = i+1;                
%                 dZdx = 1./diff(UU(2:end-1,end-1:end,1),1,2);
%                 dEdy = 1./(UU(3:end,end,2) - UU(1:end-2,end,2));
            end
            dZdx = nZeta;
            dEdy = nEta;
            
            UU(2:end-1, ind,3) = (yEta.*DxZ - xEta.*DyZ)./jac - cBC.*(dZdx);
            UU(2:end-1, ind,4) = (xZeta.*DyZ - yZeta.*DxZ)./jac - cBC.*(dEdy);
            UU(end,i,3) = (yEta(end)*DxZ(end) - xEta(end)*DyZ(end))./jac(end)...
                            - cBC.*(nZeta.*(2-4*UU(end,i,1)) + 0.5*nZeta*sqrt(2*UU(end,i,2)-1))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            if (UU(end,i,1) - UU(end-1,i,1)) ~= 0
                UU(end,i,4) = (xZeta(end)*DyZ(end) - yZeta(end)*DxZ(end))./jac(end)...
                            - cBC.*((2-4*UU(end,i,1))./(UU(end,i,1) - UU(end-1,i,1)) + 1./(UU(end,i,2)-UU(end-1,i,2)))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            else
                UU(end,i,4) = (xZeta(end)*DyZ(end) - yZeta(end)*DxZ(end))./jac(end)...
                            - cBC.*(1./(UU(end,i,2)-UU(end-1,i,2)))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            end
            UU(1,i,3) = (yEta(1)*DxZ(1) - xEta(1)*DyZ(1))./jac(1);
            UU(1,i,4) = (xZeta(1)*DyZ(1) - yZeta(1)*DxZ(1))./jac(1) - cBC.*(1./(UU(2,i,2) - UU(1,i,2)));
        else
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
                            - cBC.*(nZeta.*(2-4*UU(end,i,1)) + 0.5*nZeta*sqrt(2*UU(end,i,2)-1))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            if (UU(end,i,1) - UU(end-1,i,1)) ~= 0
                UU(end,i,4) = (xZeta(end)*DyN - yZeta(end)*DxN)./jac(end)...
                                - cBC.*((2-4*UU(end,i,1))./(UU(end,i,1) - UU(end-1,i,1)) + 1./(UU(end,i,2)-UU(end-1,i,2)))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            else
                UU(end,i,4) = (xZeta(end)*DyN - yZeta(end)*DxN)./jac(end)...
                                - cBC.*(1./(UU(end,i,2)-UU(end-1,i,2)))./(16*UU(end,i,1)^2-16*UU(end,i,1)+5);
            end
            UU(1,i,3) = (yEta(1)*DxS - xEta(1)*DyS)./jac(1);
            UU(1,i,4) = (xZeta(1)*DyS - yZeta(1)*DxS)./jac(1) - cBC.*(1./(UU(2,i,2) - UU(1,i,2)));
        end        
        
        % Calculate thomas coefficients
        aa = [repmat([0; gamma - (jac.^2).*UU(2:end-1,i,4); 0],2,1);...
                repmat([0; gamma; 0],2,1)];
        bb = repmat([1; -2.*(alpha + gamma); 1], size(UU,3),1);
        cc = [repmat([0; gamma + (jac.^2).*UU(2:end-1,i,4); 0], 2, 1);...
                repmat([0; gamma; 0],2,1)];
        dd = [  UU(1,i,1); -alpha.*sum(UU(2:end-1,[i+1,i-1],1), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],1), 1, 2)-diff(UU(1:end-2,[i+1,i-1],1), 1, 2)) - (jac.^2).*UU(2:end-1,i,3).*xZeta; UU(end,i,1);...
                UU(1,i,2); -alpha.*sum(UU(2:end-1,[i+1,i-1],2), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],2), 1, 2)-diff(UU(1:end-2,[i+1,i-1],2), 1, 2)) - (jac.^2).*UU(2:end-1,i,3).*yZeta; UU(end,i,2);...
                UU(1,i,3); -alpha.*sum(UU(2:end-1,[i+1,i-1],3), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],3), 1, 2)-diff(UU(1:end-2,[i+1,i-1],3), 1, 2)); UU(end,i,3);...
                UU(1,i,4); -alpha.*sum(UU(2:end-1,[i+1,i-1],4), 2) + 2.*beta.*(diff(UU(3:end,[i+1,i-1],4), 1, 2)-diff(UU(1:end-2,[i+1,i-1],4), 1, 2)); UU(end,i,4)];
            
        UU(:,i,:) = omega.*reshape(thomas3(aa, bb, cc, dd),size(UU,1),1,size(UU,3)) + (1-omega).*UU(:,i,:);
        
    end
    
    %% Calculate Residual
    res(end+1) = meshResidual(UU);
end