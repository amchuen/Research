function varargout = fluxCD_jst(eqnFunc, GR, FL, BC, EE)
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
[FF, GG, ~, PP1, PP2] = eqnFunc(GR, FL, BC, EE);

% Sound Speed Using Ideal-Gas Relations
V_C1 = sqrt(FL.gam.*PP1./EE1(:,:,1));
V_C2 = sqrt(FL.gam.*PP2./EE2(:,:,1));
V2 = EE2(:,:,indP2)./EE2(:,:,indRho);%./GR.RR;
V1 = EE1(:,:,indP1)./EE1(:,:,indRho);%./GR.RR;

EE2(:,:,indE) = EE2(:,:,indE) + PP2;
EE1(:,:,indE) = EE1(:,:,indE) + PP1;

%% Compute Viscosities
k2 = 1;
k4 = k2/32;
c4 = 2*k2;
% Plim = 1e-7;

r_1 = abs(V1) + V_C1;
r_N = max(r_1(3:end, :) , r_1(2:end-1,:));
r_S = max(r_1(1:end-2, :) , r_1(2:end-1,:));
r_2 = abs(V2) + V_C2;
r_E = max(r_2(:,3:end) , r_2(:,2:end-1));
r_W = max(r_2(:,1:end-2) , r_2(:,2:end-1));

% rr_1 = V1 + V_C1;
% % rr_N = 0.5.*(rr_1(3:end, :) + rr_1(2:end-1,:));
% % rr_S = 0.5.*(rr_1(1:end-2, :) + rr_1(2:end-1,:));
% rr_2 = V2 + V_C2;
% rr_E = 0.5.*(rr_2(:,3:end) + rr_2(:,2:end-1));
% rr_W = 0.5.*(rr_2(:,1:end-2) + rr_2(:,2:end-1));

% s_NS = [zeros(size(PP1(1,:))); abs((PP1(3:end,:)-2.*PP1(2:end-1,:)+PP1(1:end-2,:)))./max(abs(PP1(3:end,:) - PP1(2:end-1,:)) + abs(PP1(2:end-1,:) - PP1(1:end-2,:)),Plim); zeros(size(PP1(end,:)))];
% s_EW = [zeros(size(PP2(:,1))), abs((PP2(:,3:end)-2.*PP2(:,2:end-1)+PP2(:,1:end-2)))./max(abs(PP2(:,3:end) - PP2(:,2:end-1)) + abs(PP2(:,2:end-1) - PP2(:,1:end-2)),Plim), zeros(size(PP2(:,end)))];
s_NS = [zeros(size(PP1(1,:))); abs((PP1(3:end,:)-2.*PP1(2:end-1,:)+PP1(1:end-2,:)))./(abs(PP1(3:end,:)) + 2.*abs(PP1(2:end-1,:)) + abs(PP1(1:end-2,:))); zeros(size(PP1(end,:)))];
s_EW = [zeros(size(PP2(:,1))), abs((PP2(:,3:end)-2.*PP2(:,2:end-1)+PP2(:,1:end-2)))./(abs(PP2(:,3:end)) + 2.*abs(PP2(:,2:end-1)) + abs(PP2(:,1:end-2))), zeros(size(PP2(:,end)))];
s_N = max(s_NS(3:end,:),s_NS(2:end-1,:));
s_S = max(s_NS(1:end-2,:),s_NS(2:end-1,:));
s_E = max(s_EW(:,3:end), s_EW(:,2:end-1));
s_W = max(s_EW(:,1:end-2), s_EW(:,2:end-1));

% eps2_NS = k2.*r_1.*s_NS;
% eps2_N = 0.5.*(eps2_NS(3:end,:) + eps2_NS(2:end-1,:));
% eps2_S = 0.5.*(eps2_NS(1:end-2,:) + eps2_NS(2:end-1,:));
% eps2_EW = k2.*r_2.*s_EW;
% eps2_E = 0.5.*(eps2_EW(:,3:end) + eps2_EW(:,2:end-1));
% eps2_W = 0.5.*(eps2_EW(:,1:end-2) + eps2_EW(:,2:end-1));
eps2_N = k2.*r_N.*s_N;
eps2_S = k2.*r_S.*s_S;
eps2_E = k2.*r_E.*s_E;
eps2_W = k2.*r_W.*s_W;

eps4_N = max(0, k4.*r_N-c4.*eps2_N);
eps4_S = max(0, k4.*r_S-c4.*eps2_S);
eps4_E = max(0, k4.*r_E-c4.*eps2_E);
eps4_W = max(0, k4.*r_W-c4.*eps2_W);

%% Perform Second-Order Central Difference Operations

if GR.isPolar
    error('Work in Progress!\n');
else
    F_U = (FF(:,3:end,:) - FF(:,1:end-2,:))./(2.*GR.dx) + (GG(3:end,:,:) - GG(1:end-2,:,:))./(2.*GR.dy);
end

%% 2nd Order Dissipation

% Second-Order Diffusion Terms
d2_N = (EE1(3:end,:,:) - EE1(2:end-1,:,:));
d2_S = (EE1(2:end-1,:,:) - EE1(1:end-2,:,:));
d2_E = (EE2(:,3:end,:) - EE2(:,2:end-1,:));
d2_W = (EE2(:,2:end-1,:) - EE2(:,1:end-2,:));

% Second Derivatives
d2Ud1 = [zeros(size(EE1(1,:,:))); d2_N - d2_S; zeros(size(EE1(end,:,:)))];
d2Ud2 = [zeros(size(EE2(:,1,:))), d2_E - d2_W, zeros(size(EE2(:,end,:)))];

% Fourth-Order Diffusion Terms
d4_N = d2Ud1(3:end,:,:) - d2Ud1(2:end-1,:,:);
d4_S = d2Ud1(2:end-1,:,:) - d2Ud1(1:end-2,:,:);
d4_E = d2Ud2(:,3:end,:) - d2Ud2(:,2:end-1,:);
d4_W = d2Ud2(:,2:end-1,:) - d2Ud2(:,1:end-2,:);

% Total Dissipation
D_U2 = (eps2_N.*d2_N - eps2_S.*d2_S)./GR.dy + (eps2_E.*d2_E - eps2_W.*d2_W)./GR.dx;
D_U4 = (eps4_N.*d4_N - eps4_S.*d4_S)./GR.dy + (eps4_E.*d4_E - eps4_W.*d4_W)./GR.dx;

D_U = D_U2 - D_U4;

flux = F_U - D_U;

%% Outputs

waveSpd = [max(abs(r_2(:))), max(abs(r_1(:)))];

varargout{1} = flux;
varargout{2} = waveSpd;
varargout{3} = V_C1(2:end-1,:);
varargout{4} = 0.5.*(eps2_N+eps2_S);
varargout{5} = 0.5.*(eps2_W+eps2_E);

end