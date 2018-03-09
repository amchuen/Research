function [vrDiff, timeCoeffs] = VRdiffusion(EE, GR, BC, FL)
% indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(UU.fv,3)); 
indV1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(EE,3));
indV2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(EE,3));

% Compute Viscosities
[epsN, epsS, epsE, epsW] = vonNeumRichtVisc(EE, GR, BC, FL);

% Compute Diffusive Terms
f2 = epsE.*[EE(:,2:end,:),bcCalc(GR, BC, EE(:,:,:), 'E')]  + epsW.*[bcCalc(GR, BC, EE(:,:,:), 'W'), EE(:,1:end-1,:)];
if GR.isPolar
    f1(2:end-1,:,:) = GR.RR_N(2:end-1,:).*EE(3:end,:,:) + GR.RR_S(2:end-1,:).*EE(1:end-2,:,:);
    f1(1,:,:) = GR.RR_N(1,:).*EE(2,:,:) + GR.RR_S(1,:).*bcCalc(GR, BC, EE(:,:,:),'S');
    f1(end,:,:) = GR.RR_S(end,:).*EE(end-1,:,:) + GR.RR_N(end,:).*bcCalc(GR,BC, EE(:,:,:),'N');

    % Calculate Vector Laplacian Terms
    if any(indV1) || any(indV2)
        rotLaplace =- indV1.*(EE(:,:,indV1) +...
                            2.*([EE(:,2:end,indV2),bcCalc(GR,BC,EE(:,:,:),'E',find(indV2))] - [bcCalc(GR,BC,EE(:,:,:),'W',find(indV2)),EE(:,1:end-1,indV2)])./(2.*GR.dT))./(GR.RR.^2)...
                - indV2.*(EE(:,:,indV2) -...
                            2.*([EE(:,2:end,indV1),bcCalc(GR,BC,EE(:,:,:),'E',find(indV1))] - [bcCalc(GR,BC,EE(:,:,:),'W',find(indV1)),EE(:,1:end-1,indV1)])./(2.*GR.dT))./(GR.RR.^2);
    else
        rotLaplace = 0;
    end
    alpha2 = 2.*GR.dt.*epsFunc(GR,BC,'T')./(GR.RR.*GR.dT).^2;
    alpha1 = 2.*GR.dt.*epsFunc(GR,BC,'R')./(GR.RR.*GR.dR.^2).*0.5.*(GR.RR_N+GR.RR_S);
else
    % Compute diffusion terms at time n
    f1 = epsN.*[EE(2:end,:,:); bcCalc(GR,BC, EE(:,:,:),'N')] + epsS.*[bcCalc(GR, BC, EE(:,:,:),'S'); EE(1:end-1,:,:)];
    vrDiff = f2./(GR.dx.^2)+f1./(GR.dy.^2);
    
    % time-stepping coefficients -> Dufort Frankel Formulation
    alpha2 = GR.dt.*(epsW + epsE)./(GR.dx.^2);
    alpha1 = GR.dt.*(epsN + epsS)./(GR.dy.^2);
end

% time-stepping coefficients -> Dufort Frankel Formulation
alphaTot = alpha2 + alpha1;
timeCoeffs{1} = 1 + 0.5.*(1+1/GR.cflFactor).*alphaTot;
timeCoeffs{2} = 0;
timeCoeffs{3} = 1 - 0.5.*(1+1/GR.cflFactor).*alphaTot;


end