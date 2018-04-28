function varargout = fluxCD_24Diff(eqnFunc, viscFunc, GR, FL, BC, EE)
% This function computes the spatial deriatives (flux) in the Euler
% equations using central differences. The inputs consist of the function
% to compute the nonlinear terms, as well as grid information. Boundary
% information is not needed as it is assumed to be computed/accounted for
% in the nonlinear function.

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(EE,3));
indE = reshape(strcmp(BC.N.varName, '\rho e'),1,1,size(EE,3));

EE2 = [bcCalc(GR,FL,BC,EE,'W'), EE, bcCalc(GR,FL,BC,EE,'E')];
EE1 = [bcCalc(GR,FL,BC,EE,'S'); EE; bcCalc(GR,FL,BC,EE,'N')];

%% Compute Flux Terms

% Nonlinear Terms from Conservation Equations
% [FF, GG, PP] = eqnFunc(GR, FL, BC, EE);
[FF, GG, ~, PP1, PP2] = eqnFunc(GR, FL, BC, EE);

% Characteristic Speeds (for the highest ones?)
C_speed2 = abs(EE2(:,:,indP2)./EE2(:,:,indRho)) + sqrt(FL.gam.*PP2./EE2(:,:,1));
C_speed1 = abs(EE1(:,:,indP1)./EE1(:,:,indRho)) + sqrt(FL.gam.*PP1./EE1(:,:,1));

% Enforce CFL
waveSpd = [max(abs(C_speed2(:))), max(abs(C_speed1(:)))];

%% Perform Second-Order Central Difference Operations

if GR.isPolar
    error('Work in Progress!\n');
else
    F_U = (FF(:,3:end,:) - FF(:,1:end-2,:))./(2.*GR.dx) + (GG(3:end,:,:) - GG(1:end-2,:,:))./(2.*GR.dy);
end

%% 2nd Order Dissipation

% Diffusion Terms
d2_N = (EE1(3:end,:,:) - EE1(2:end-1,:,:));%./GR.dy;
d2_S = (EE1(2:end-1,:,:) - EE1(1:end-2,:,:));%./GR.dy;
d2_E = (EE2(:,3:end,:) - EE2(:,2:end-1,:));%./GR.dx;
d2_W = (EE2(:,2:end-1,:) - EE2(:,1:end-2,:));%./GR.dx;

% Compute Viscosities
[eps2_N, eps2_S, eps2_E, eps2_W] = viscFunc(EE, GR, BC, FL);

% Total Dissipation
D_U2 = (eps2_N.*d2_N - eps2_S.*d2_S)./(GR.dy^2) + (eps2_E.*d2_E - eps2_W.*d2_W)./(GR.dx^2);

%% Fourth Order Dissipation
% Second Derivatives
d2Ud1 = [zeros(size(EE1(1,:,:))); d2_N - d2_S; zeros(size(EE1(end,:,:)))];
d2Ud2 = [zeros(size(EE2(:,1,:))), d2_E - d2_W, zeros(size(EE2(:,end,:)))];

% Fourth-Order Diffusion Terms
% d4_N = d2Ud1(3:end,:,:) - d2Ud1(2:end-1,:,:);
% d4_S = d2Ud1(2:end-1,:,:) - d2Ud1(1:end-2,:,:);
% d4_E = d2Ud2(:,3:end,:) - d2Ud2(:,2:end-1,:);
% d4_W = d2Ud2(:,2:end-1,:) - d2Ud2(:,1:end-2,:);

% Compute Viscosities
incFact = 0.0625;
eps4_N = max(0,(incFact.*max(C_speed1(3:end,:),C_speed1(2:end-1,:))./8 - 0.25.*eps2_N./(GR.dy)));
eps4_S = max(0,(incFact.*max(C_speed1(1:end-2,:),C_speed1(2:end-1,:))./8 - 0.25.*eps2_S./(GR.dy)));
eps4_E = max(0,(incFact.*max(C_speed2(:,3:end),C_speed2(:,2:end-1))./8 - 0.25.*eps2_E./(GR.dx)));
eps4_W = max(0,(incFact.*max(C_speed2(:,1:end-2),C_speed2(:,2:end-1))./8 - 0.25.*eps2_W./(GR.dx)));

% Compute Dissipation
D_U4 = (eps4_N.*(d2Ud1(3:end,:,:) - d2Ud1(2:end-1,:,:)) - eps4_S.*(d2Ud1(2:end-1,:,:) - d2Ud1(1:end-2,:,:)))./GR.dy + (eps4_E.*(d2Ud2(:,3:end,:) - d2Ud2(:,2:end-1,:)) - eps4_W.*(d2Ud2(:,2:end-1,:) - d2Ud2(:,1:end-2,:)))./GR.dx;


%% Calculate Total Flux with Dissipation

flux = F_U - D_U2 + D_U4;

%% Outputs
varargout{1} = flux;
varargout{2} = waveSpd;
varargout{3} = sqrt(FL.gam.*PP2(:,2:end-1)./EE2(:,2:end-1,1));
varargout{4} = 0.5.*(eps2_N+eps2_S);
varargout{5} = 0.5.*(eps2_W+eps2_E);

end