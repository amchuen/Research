function [fvOUT, waveSpd] = eulerIsenFunc(GR, FL, BC, EE)

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varType, 's'),1,1,size(EE,3));

% Pressure Calculation
PP = EE(:,:,indRho).^(FL.gam)./(FL.gam .* FL.M0.^2);

% Boundary Calculation -> include pressure terms
bcN = bcCalc(GR,BC,EE,'N') .* bcCalc(GR,BC,EE,'N',indP1) ./ bcCalc(GR,BC,EE,'N',indRho);bcNP = indP1.*(bcCalc(GR,BC,EE,'N',indRho).^(FL.gam))./(FL.gam.*FL.M0.^2);
bcS = bcCalc(GR,BC,EE,'S') .* bcCalc(GR,BC,EE,'S',indP1) ./ bcCalc(GR,BC,EE,'S',indRho);bcSP = indP1.*(bcCalc(GR,BC,EE,'S',indRho).^(FL.gam))./(FL.gam.*FL.M0.^2);
bcE = bcCalc(GR,BC,EE,'E') .* bcCalc(GR,BC,EE,'E',indP2) ./ bcCalc(GR,BC,EE,'E',indRho);bcEP = indP2.*(bcCalc(GR,BC,EE,'E',indRho).^(FL.gam))./(FL.gam.*FL.M0.^2);
bcW = bcCalc(GR,BC,EE,'W') .* bcCalc(GR,BC,EE,'W',indP2) ./ bcCalc(GR,BC,EE,'W',indRho);bcWP = indP2.*(bcCalc(GR,BC,EE,'W',indRho).^(FL.gam))./(FL.gam.*FL.M0.^2);

% Calculate Wave Terms (Conservation Form)
FF = EE.*EE(:,:,indP2)./EE(:,:,indRho);% + indP2.*PP;
GG = EE.*EE(:,:,indP1)./EE(:,:,indRho);% + indP1.*PP;

if GR.isPolar
    %FF = FF + indP2.*PP;
    % radial derivatives must calculate pressure and vel. separately!
    GG_1 =  (GR.RR_N.*[GG(2:end,:,:);bcN] - GR.RR_S.*[bcS;GG(1:end-1,:,:)])./(2.*GR.dR.*GR.RR); ... % vector Flux
    PP_1 =  ([indP1.*PP(2:end,:);bcNP] - [bcSP;indP1.*PP(1:end-1,:)])./(2*GR.dR);...Pressure Term
    C_1 =  -indP1.*(EE(:,:,indP2).^2)./(GR.RR.*EE(:,:,indRho)); ... Centrifugal Term
%     FF_2 =  ([FF(:,2:end,:),bcE+bcEP] - [bcW+bcWP,FF(:,1:end-1,:)])./(2.*GR.dT.*GR.RR)...
    FF_2 =  ([FF(:,2:end,:),bcE] - [bcW,FF(:,1:end-1,:)])./(2.*GR.dT.*GR.RR);...
    PP_2 = ([indP2.*PP(:,2:end),bcEP] - [bcWP,indP2.*PP(:,1:end-1)])./(2.*GR.dT.*GR.RR);
    C_2 = indP2.*EE(:,:,indP1).*EE(:,:,indP2)./(EE(:,:,indRho).*GR.RR);

    fvOUT = FF_2 + GG_1 + PP_1 + PP_2 + C_1 + C_2;
else
    bcN = bcN + bcNP; bcS = bcS + bcSP;
    bcE = bcE + bcEP; bcW = bcW + bcWP;
    FF = FF + indP2.*PP;
    GG = GG + indP1.*PP;
    FF_2 = ([FF(:,2:end,:),bcE] - [bcW,FF(:,1:end-1,:)])./(2.*GR.dx);
    GG_1 = ([GG(2:end,:,:);bcN] - [bcS;GG(1:end-1,:,:)])./(2.*GR.dy);
    fvOUT = FF_2 + GG_1;
end

% rho = FF(:,:,1);
phi2w = EE(:,:,indP2)./EE(:,:,indRho);%./GR.RR;
phi1w = EE(:,:,indP1)./EE(:,:,indRho);%./GR.RR;
% rhow = rho;%./GR.RR;
% waveSpd = [max([max(abs(phi2w(:))), max(rhow(:))]), max([max(abs(phi1w(:))), max(rhow(:))])];
waveSpd = [max(abs(phi2w(:))), max(abs(phi1w(:)))];

end