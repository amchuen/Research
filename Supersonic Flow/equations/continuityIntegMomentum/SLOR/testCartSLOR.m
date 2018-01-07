clc;
close all;
clear;

%% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
% Define Grid
dx = 0.05;
dy = 0.08;

% Field Axis Values
y_max = 50;%.*dy; %25;
% x_min = -20*dx;
x_max = 11;
% x_max = 39*dx; %(-19*dx);
x_min = -10;
x_vals = x_min:dx:x_max;
y_vals = 0:dy:y_max;

[XX, YY] = meshgrid(x_vals, y_vals);

GR.dx = dx;
GR.dy = dy;
GR.XX = XX;
GR.YY = YY;

%% Define Fluid Parameters
gam = 1.4;
M0 = 0.85;

%% Simulation Parameters
tol = 1e-5;
w = 1;

%% Diffusion Coefficients

epsFunc = @(GR, BC, DIR) 0.005;

%% Boundary Conditions

% % Body Values - Ramp
tau = 0.05;
m_x = tand(8); % dy/dx
% x_vals = x_vals;
% dx = dx;
YY_B = [zeros(size(x_vals(x_vals <0))), ...
        2*tau.*x_vals((x_vals>=0)&(x_vals <=1)).*(1- x_vals((x_vals>=0)&(x_vals <=1))),...
        zeros(size(x_vals(x_vals >1)))];
dyBdx = zeros(size(YY_B));

for i = 2:(length(YY_B)-1)
   dyBdx(i) = (YY_B(i+1) - YY_B(i-1))/(2*dx);
end

% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)

% Far-field... need to update this?
BC.N.physical = 'inlet';
BC.N.val = {1, x_vals};
BC.N.varType = {'s','phi'};
BC.N.dydx = 0;

% Inlet
BC.W = BC.N;
BC.W.val = {1, x_vals(1)};

% Wall
BC.S.physical = 'wall';
BC.S.varType = BC.N.varType;
BC.S.dydx = dyBdx;

% Outlet
% BC.E.physical = 'outlet';
% BC.E.varType = BC.N.varType;
BC.E = BC.N;
BC.E.val = {1, x_vals(end)};

%% Run Simulation

FV = cat(3, ones(size(GR.XX)), GR.XX);

res1 = [1, 1];
%while norm(res1(end,:)) > tol
for ii = 1:500
    FV(:,1,:) = bcCalc(GR,BC,FV,'W');
    FV(:,end,:) = bcCalc(GR,BC,FV,'E');
    resiter = zeros(2,size(FV,2)-2);
	
    for i = 2:size(FV,2)-1
        res2 = [1,1];
        while norm(res2) > tol
            % Apply boundary conditions
            FV(end,:,:) = bcCalc(GR,BC,FV, 'N');
            FV(1,:,:) = bcCalc(GR,BC,FV, 'S');

            % Calculate velocities
            phiy_n = diff(FV(2:end,i,1))./dy;
            phiy_s = diff(FV(1:end-1,i,1))./dy;
            phix_w = (FV(2:end-1,i,1)-FV(2:end-1,i-1,1))./dx;
            phix_e = (FV(2:end-1,i+1,1)-FV(2:end-1,i,1))./dx;

            % Set up equation systems
            aRho = [0; epsFunc(GR,BC,'Y')./dy^2 + phiy_s./(2*dy); 0];
            bRho = [1; -(2*(epsFunc(GR,BC,'X')./(dx^2) + epsFunc(GR,BC,'Y')./(dy^2))+ 0.5.*(phiy_n-phiy_s)./dy + 0.5.*(phix_e - phix_w)./dx);1];
            cRho = [0; epsFunc(GR,BC,'Y')./dy^2 - phiy_n./(2*dy); 0];
            dRho = [FV(1,i,1);  (-epsFunc(GR,BC,'X')./(dx^2)+phix_e./(2*dx)).*FV(2:end-1,i+1,1)...
                                + (-epsFunc(GR,BC,'X')./(dx^2)-phix_w./(2*dx)).*FV(2:end-1,i-1,1); FV(end,i,1)];

            aPhi = [0; epsFunc(GR,BC,'Y')./(dy^2).*ones(size(FV,1)-2,1); 0];
            bPhi = [1; (-2.*(epsFunc(GR,BC,'X')./(dx^2) + epsFunc(GR,BC,'Y')./(dy^2))).*ones(size(FV(2:end-1,i,2))); 1];
            cPhi = aPhi;
            dPhi = [FV(1,i,2);  - epsFunc(GR,BC,'X')./(dx^2).*(FV(2:end-1,i+1,2)+FV(2:end-1,i-1,2))...
                                + 0.5.*((0.5.*(phiy_n+phiy_s)).^2+(0.5.*(phix_w+phix_e)).^2 - 1)...
                                + (FV(2:end-1,i,1).^(gam-1)-1)./((gam-1).*M0.^2); FV(end,i,2)];

            fvStar = thomas3([aRho;aPhi], [bRho;bPhi], [cRho;cPhi], [dRho; dPhi]);
            corr = cat(3,fvStar(1:size(FV,1)), fvStar(size(FV,1)+1:end)) - FV(:,i,:);
            FV(:,i,:) = FV(:,i,:) + w.*corr;
            
            resRho = abs(aRho(2:end-1).*FV(1:end-2,i,1)+bRho(2:end-1).*FV(2:end-1,i,1)+cRho(2:end-1).*FV(3:end,i,1) - dRho(2:end-1));
            resPhi = abs(aPhi(2:end-1).*FV(1:end-2,i,2)+bPhi(2:end-1).*FV(2:end-1,i,2)+cPhi(2:end-1).*FV(3:end,i,2) - dPhi(2:end-1));
            res2 = [max(abs(resRho(:))), max(abs(resPhi))];
        end
        resiter(:, i-1) = res2';

    end
        
    res1(end+1,:) = [max(resiter(1,:)), max(resiter(2,:))];

    figure(1);
%     semilogy(1:length(res), res(:,1));
%     hold on; semilogy(1:length(res), res(:,2));
%     legend('\rho', '\phi', 'Location', 'Best');
    contourf(XX, YY, FV(:,:,1),50);
%     title('Residual');
%     hold off;
%     figure(2);
%     contourf(XX,YY,FV(:,:,2),50);
    %colorbar;
    drawnow;

end