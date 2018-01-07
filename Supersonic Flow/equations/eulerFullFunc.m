function [fvOUT, waveSpd] = eulerFullFunc(GR, FL, BC, EE)

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(EE,3));
indE = reshape(strcmp(BC.N.varName, '\rho e'),1,1,size(EE,3));

% Pressure Calculation
PP = EE(:,:,indRho).^(FL.gam)./(FL.gam .* FL.M0.^2);

% Boundary Calculation -> include pressure terms
bcN = (bcCalc(GR,BC,EE,'N') + indE.*bcCalc(GR,BC,EE,'N',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2)).* bcCalc(GR,BC,EE,'N',indP1) ./ bcCalc(GR,BC,EE,'N',indRho) + indP1.*bcCalc(GR,BC,EE,'N',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2);
bcS = (bcCalc(GR,BC,EE,'S') + indE.*bcCalc(GR,BC,EE,'S',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2)).* bcCalc(GR,BC,EE,'S',indP1) ./ bcCalc(GR,BC,EE,'S',indRho) + indP1.*bcCalc(GR,BC,EE,'S',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2);
bcE = (bcCalc(GR,BC,EE,'E') + indE.*bcCalc(GR,BC,EE,'E',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2)).* bcCalc(GR,BC,EE,'E',indP2) ./ bcCalc(GR,BC,EE,'E',indRho) + indP2.*bcCalc(GR,BC,EE,'E',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2);
bcW = (bcCalc(GR,BC,EE,'W') + indE.*bcCalc(GR,BC,EE,'W',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2)).* bcCalc(GR,BC,EE,'W',indP2) ./ bcCalc(GR,BC,EE,'W',indRho) + indP2.*bcCalc(GR,BC,EE,'W',indRho).^(FL.gam)./(FL.gam.*FL.M0.^2);

% Calculate Wave Terms (Conservation Form)
FF = (EE+indE.*PP).*EE(:,:,indP2)./EE(:,:,indRho) + indP2.*PP;
GG = (EE+indE.*PP).*EE(:,:,indP1)./EE(:,:,indRho) + indP1.*PP;

FF_2 = ([FF(:,2:end,:),bcE] - [bcW,FF(:,1:end-1,:)])./(2.*GR.dx);
GG_1 = ([GG(2:end,:,:);bcN] - [bcS;GG(1:end-1,:,:)])./(2.*GR.dy);

fvOUT = FF_2 + GG_1;

% rho = FF(:,:,1);
phi2w = EE(:,:,indP2)./EE(:,:,indRho);%./GR.RR;
phi1w = EE(:,:,indP1)./EE(:,:,indRho);%./GR.RR;
% rhow = rho;%./GR.RR;
% waveSpd = [max([max(abs(phi2w(:))), max(rhow(:))]), max([max(abs(phi1w(:))), max(rhow(:))])];
waveSpd = [max(abs(phi2w(:))), max(abs(phi1w(:)))];

end