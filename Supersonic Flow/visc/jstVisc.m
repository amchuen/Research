function [jstDiff, timeCoeffs] = jstVisc(FF, GR, BC, FL, varargin)

%% Calculate Viscosities

% Coefficients
K_4 = 1/256;
K_2 = 0.25;

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(FF,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(FF,3));
indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(FF,3));

% Calculate Pressure Sensor
PP_WE = ([bcCalc(GR,BC,FF,'W',indRho), FF, bcCalc(GR,BC,FF,'E', indRho)].^FL.gam)./(FL.gam.*FL.M0^2);
PP_NS = ([bcCalc(GR,BC,FF,'S', indRho); FF; bcCalc(GR,BC,FF,'N', indRho)].^FL.gam)./(FL.gam.*FL.M0^2);
nu_WE = abs(PP_WE(:,3:end) - 2*PP_WE(:,2:end-1) + PP_WE(:,1:end-2))./(abs(PP_WE(:,3:end)) + 2.*abs(PP_WE(:,2:end-1)) + abs(PP_WE(:,1:end-2)));
nu_NS = abs(PP_NS(3:end,:) - 2*PP_NS(:,2:end-1) + PP_NS(1:end-2,:))./(abs(PP_NS(3:end,:)) + 2.*abs(PP_NS(2:end-1,:)) + abs(PP_NS(1:end-2,:)));

% Repeat Terms at Boundaries?
nu_WE = [nu_WE(:,1), nu_WE, nu_WE(:,end)];
nu_NS = [nu_NS(1,:); nu_NS, nu_NS(end,:)];

% Calculate epsilon-2's
eps2_E = K_2.*max(nu_WE(:,3:end),nu_WE(:,2:end-1));
eps2_W = K_2.*max(nu_WE(:,1:end-2),nu_WE(:,2:end-1));
eps2_N = K_2.*max(nu_NS(3:end,:),nu_NS(2:end-1,:));
eps2_S = K_2.*max(nu_NS(1:end-2,:),nu_NS(2:end-1,:));

% Calculate epsilon-4's
eps4_E = max(0, (K_4 - eps2_E));
eps4_W = max(0, (K_4 - eps2_W));
eps4_N = max(0, (K_4 - eps2_N));
eps4_S = max(0, (K_4 - eps2_S));

%% Calculate 2nd and 4th Differences 

% 2nd
diff2 = eps2_E.*[FF(:,2:end,:),bcCalc(GR, BC, FF(:,:,:), 'E')]  + eps2_W.*[bcCalc(GR, BC, FF(:,:,:), 'W'), FF(:,1:end-1,:)] + ...
        eps2_N.*[FF(2:end,:,:); bcCalc(GR,BC, FF(:,:,:),'N')] + eps2_S.*[bcCalc(GR, BC, FF(:,:,:),'S'); FF(1:end-1,:,:)];

% 4th

end