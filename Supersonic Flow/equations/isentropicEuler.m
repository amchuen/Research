function [FF, GG, PP, varargout] = isentropicEuler(GR, FL, BC, EE)

% Get Indexing
indP1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indP2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));
indRho = reshape(strcmp(BC.N.varType, 's'),1,1,size(EE,3));

% Pressure Calculation
PP = (EE(:,:,indRho).^FL.gam)./(FL.gam .* FL.M0.^2);

% Nonlinear Transport Terms
EEx = [bcCalc(GR,FL,BC,EE,'W'), EE, bcCalc(GR,FL,BC,EE,'E')];
EEy = [bcCalc(GR,FL,BC,EE,'S'); EE; bcCalc(GR,FL,BC,EE,'N')];

FF =    indRho.*EEx(:,:,indP2)...
        + indP2.*(EEx(:,:,indP2).^2./EEx(:,:,indRho) + (EEx(:,:,indRho).^FL.gam)./(FL.gam.*FL.M0^2))...
        + indP1.*(EEx(:,:,indP1).*EEx(:,:,indP2)./EEx(:,:,indRho));
GG =    indRho.*EEy(:,:,indP1)...
        + indP1.*(EEy(:,:,indP1).^2./EEy(:,:,indRho) + (EEy(:,:,indRho).^FL.gam)./(FL.gam.*FL.M0^2))...
        + indP2.*(EEy(:,:,indP1).*EEy(:,:,indP2)./EEy(:,:,indRho));
    
% More Pressure Terms for JST Scheme
PP2 = (EEx(:,:,indRho).^FL.gam)./(FL.gam .* FL.M0.^2);
PP1 = (EEy(:,:,indRho).^FL.gam)./(FL.gam .* FL.M0.^2);

varargout{1} = PP1;
varargout{2} = PP2;
end