function [fvOUT, waveSpd] = centFD(GR, FL, BC, FF, eqnFunc)
% This function computes the spatial deriatives (flux) in the Euler
% equations using central differences. The inputs consist of the function
% to compute the nonlinear terms, as well as grid information. Boundary
% information is not needed as it is assumed to be computed/accounted for
% in the nonlinear function.


%% Compute Flux Terms

[fluxVals, waveSpd] = eqnFunc(GR, FL, BC, FF);

%% Perform Second-Order Central Difference Operations

if GR.isPolar
    error('Work in Progress!\n');
    
else
    fvOUT = (fluxVals.FF(:,3:end,:) - fluxVals.FF(:,1:end-2,:))./(2.*GR.dx) + ...
            (fluxVals.GG(3:end,:,:) - fluxVals.GG(1:end-2,:,:))./(2.*GR.dy);
    
end

end