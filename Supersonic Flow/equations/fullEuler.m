function [FF, GG, PP, varargout] = fullEuler(GR, FL, BC, EE)

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(EE,3));
indE = reshape(strcmp(BC.N.varName, '\rho e'),1,1,size(EE,3));

% Pressure Calculation
PP = (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho));

% Nonlinear Transport Terms
EEx = [bcCalc(GR,FL,BC,EE,'W'), EE, bcCalc(GR,FL,BC,EE,'E')];
EEy = [bcCalc(GR,FL,BC,EE,'S'); EE; bcCalc(GR,FL,BC,EE,'N')];

FF =    indRho.*EEx(:,:,indP2)...
        + indP2.*(EEx(:,:,indP2).^2./EEx(:,:,indRho) + (FL.gam-1).*(EEx(:,:,indE)-0.5.*(EEx(:,:,indP1).^2 + EEx(:,:,indP2).^2)./EEx(:,:,indRho)))...
        + indP1.*(EEx(:,:,indP1).*EEx(:,:,indP2)./EEx(:,:,indRho))...
        + indE.*(EEx(:,:,indP2)./EEx(:,:,indRho)).*(EEx(:,:,indE) + (FL.gam-1).*(EEx(:,:,indE)-0.5.*(EEx(:,:,indP1).^2 + EEx(:,:,indP2).^2)./EEx(:,:,indRho)));
    
GG =    indRho.*EEy(:,:,indP1)...
        + indP1.*(EEy(:,:,indP1).^2./EEy(:,:,indRho) + (FL.gam-1).*(EEy(:,:,indE)-0.5.*(EEy(:,:,indP1).^2 + EEy(:,:,indP2).^2)./EEy(:,:,indRho)))...
        + indP2.*(EEy(:,:,indP1).*EEy(:,:,indP2)./EEy(:,:,indRho))...
        + indE.*(EEy(:,:,indP1)./EEy(:,:,indRho)).*(EEy(:,:,indE) + (FL.gam-1).*(EEy(:,:,indE)-0.5.*(EEy(:,:,indP1).^2 + EEy(:,:,indP2).^2)./EEy(:,:,indRho)));
    
% More Pressure Terms for JST Scheme
PP2 = (FL.gam-1).*(EEx(:,:,indE)-0.5.*(EEx(:,:,indP1).^2 + EEx(:,:,indP2).^2)./EEx(:,:,indRho));
PP1 = (FL.gam-1).*(EEy(:,:,indE)-0.5.*(EEy(:,:,indP1).^2 + EEy(:,:,indP2).^2)./EEy(:,:,indRho));

varargout{1} = PP1;
varargout{2} = PP2;
end