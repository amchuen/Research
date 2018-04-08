function varargout = fluxCD_2Diff(eqnFunc, viscFunc, GR, FL, BC, EE)
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

%% Compute Flux Terms

% Nonlinear Terms from Conservation Equations
[FF, GG, PP] = eqnFunc(GR, FL, BC, EE);

% Sound Speed Using Ideal-Gas Relations
V_C = sqrt(FL.gam.*PP./EE(:,:,1));
V2 = EE(:,:,indP2)./EE(:,:,indRho);%./GR.RR;
V1 = EE(:,:,indP1)./EE(:,:,indRho);%./GR.RR;
V2 = max(cat(3, abs(V2+V_C), abs(V2), abs(V2-V_C)),[],3);
V1 = max(cat(3, abs(V1+V_C), abs(V1), abs(V1-V_C)),[],3);
waveSpd = [max(abs(V2(:))), max(abs(V1(:)))];

%% Perform Second-Order Central Difference Operations

if GR.isPolar
    error('Work in Progress!\n');
else
    F_U = (FF(:,3:end,:) - FF(:,1:end-2,:))./(2.*GR.dx) + (GG(3:end,:,:) - GG(1:end-2,:,:))./(2.*GR.dy);
end

%% 2nd Order Dissipation
EEx = [bcCalc(GR,FL,BC,EE,'W'), EE, bcCalc(GR,FL,BC,EE,'E')];
EEy = [bcCalc(GR,FL,BC,EE,'S'); EE; bcCalc(GR,FL,BC,EE,'N')];

% Diffusion Terms
d2_N = (EEy(3:end,:,:) - EEy(2:end-1,:,:))./GR.dy;
d2_S = (EEy(2:end-1,:,:) - EEy(1:end-2,:,:))./GR.dy;
d2_E = (EEx(:,3:end,:) - EEx(:,2:end-1,:))./GR.dx;
d2_W = (EEx(:,2:end-1,:) - EEx(:,1:end-2,:))./GR.dx;

% Compute Viscosities
[epsN, epsS, epsE, epsW] = viscFunc(EE, GR, BC, FL);

% Total Dissipation
D_U = (epsN.*d2_N - epsS.*d2_S)./GR.dy + (epsE.*d2_E - epsW.*d2_W)./GR.dx;

flux = F_U - D_U;

%% Outputs
varargout{1} = flux;
varargout{2} = waveSpd;
varargout{3} = V_C;
varargout{4} = 0.5.*(epsN+epsS);
varargout{5} = 0.5.*(epsW+epsE);

end