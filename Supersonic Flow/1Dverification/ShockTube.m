clc;
clear;
close all;

%% Generate Grid

% dx = 0.0025;
% xx = (0:dx:1)';
xx = linspace(0,1,401)';
dx = xx(2) - xx(1);
dt = 0.01;

Ldx = diag([0, ones(1,length(xx)-2)./2.*dx],1) - diag([ones(1,length(xx)-2)./2.*dx, 0], -1);
Ldx2 = diag([0, ones(1,length(xx)-2)./(dx.^2)],1) + diag([ones(1,length(xx)-2)./(dx.^2), 0], -1);
Cdx2 = diag([0, ones(1,length(xx)-2)./(dx^2),0]);

%% Fluid Properties
gam = 1.4;
visc = dx*0.007;
cfl = 0.75;

%% Initialize Field Values

UU = ones([3, length(xx)]);

% density
UU(1,:) = (xx<=0.5)+0.125.*(xx>0.5); 

% velocity
UU(2,:) = 0; 

%Energy
PP = (xx<=0.5)+0.1.*(xx>0.5);
UU(3,:) = PP'./(gam - 1) + 0.5.*UU(2,:).^2./UU(1,:);

UU = repmat(UU',1,1,3);

%% Setup Time Simulation
time = 0;
tEnd = 10*5;

%% Run Simulation

while time(end) < tEnd
    
    % Update Field Values
    UU(:,:,1:2) = UU(:,:,2:3);
    
    % Update Flux
    PP = (gam - 1).*(UU(:,3,2) - 0.5.*UU(:,2,2).^2./UU(:,1,2));
    FF =  [ UU(:,2,2),...
            UU(:,2,2).^2./UU(:,1,2)+ PP,...
            (UU(:,3,2)+PP).*UU(:,2,2)./UU(:,1,2)];
    
    % Check CFL
    epsilon = visc.* abs((Ldx2 - 2.*Cdx2)*UU(:,1,2))./((Ldx2+2.*Cdx2)*UU(:,1,2));
    epsilon = [0; epsilon(2:end-1); 0];
    alphaL = 2.*dt.*epsilon.*Ldx2;
    alphaC = 2.*dt.*epsilon.*Cdx2;
    CFL = (UU(epsilon~=0,2,2)./UU(epsilon~=0,1,2)).*dt./dx;
    CFL = max(abs(CFL(:)));
    cflFactor = 1;
    if (~isempty(CFL) && (CFL~=cfl)) && any(epsilon > 0)
        cflFactor = min(cfl/CFL,cfl/max(abs(diag(alphaC))));
        alphaL = alphaL.*(1+cflFactor)./2;
        alphaC = alphaC.*(1+cflFactor)./2;
        dt = dt.*cflFactor;
        
%         if abs(log10(cflFactor)) > 0.005
%             fprintf('CFL condition not met!\n');
%             fprintf('Changing time steps\n');
%             fprintf('New time step:%0.5e\n', dt);
%         end
        
    end
    waveSpd = (1+cflFactor).*dt.*Ldx;

    
    UU(:,:,3) = (eye(length(xx))+alphaC)\((eye(length(xx)) - alphaC)*UU(:,:,1) + alphaL*UU(:,:,2) - waveSpd*FF);
    time(end+1) = time(end)+dt;
end

%% Post Process

figure();
plot(xx,UU(:,1,end));
title('Density');

figure();
plot(xx,(UU(:,2,end)./UU(:,1,end)));
title('Velocity');

figure();
PP = (gam - 1).*(UU(:,3,3) - 0.5.*UU(:,2,3).^2./UU(:,1,3));
plot(xx,PP);
title('Pressure');