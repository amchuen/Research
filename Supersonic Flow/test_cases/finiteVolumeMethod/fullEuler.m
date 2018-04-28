function [FF, GG, PP, varargout] = fullEuler(FL, BC, EE)

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(EE,3));
indE = reshape(strcmp(BC.N.varName, '\rho e'),1,1,size(EE,3));

% Pressure Calculation
PP = (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho));

% Nonlinear Transport Terms
FF =    indRho.*EE(:,:,indP2)...
        + indP2.*(EE(:,:,indP2).^2./EE(:,:,indRho) + (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho)))...
        + indP1.*(EE(:,:,indP1).*EE(:,:,indP2)./EE(:,:,indRho))...
        + indE.*(EE(:,:,indP2)./EE(:,:,indRho)).*(EE(:,:,indE) + (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho)));
    
GG =    indRho.*EE(:,:,indP1)...
        + indP1.*(EE(:,:,indP1).^2./EE(:,:,indRho) + (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho)))...
        + indP2.*(EE(:,:,indP1).*EE(:,:,indP2)./EE(:,:,indRho))...
        + indE.*(EE(:,:,indP1)./EE(:,:,indRho)).*(EE(:,:,indE) + (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho)));
    
% More Pressure Terms for JST Scheme
% PP2 = (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho));
% PP1 = (FL.gam-1).*(EE(:,:,indE)-0.5.*(EE(:,:,indP1).^2 + EE(:,:,indP2).^2)./EE(:,:,indRho));
% 
% varargout{1} = PP1;
% varargout{2} = PP2;
end