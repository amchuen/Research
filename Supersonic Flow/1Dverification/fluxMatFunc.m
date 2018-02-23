function flux = fluxMatFunc(Uvec, gam, dx)

% Inputs:
% Uvec - vector of all fluid parameters
% gam - heat coefficient ratio

% Procedure
% 1) Initialize and store nnz's into a column

nNodes = length(Uvec)/3;
varInd = (1:nNodes:length(Uvec))-1; % starting index for each eqn
foreInd = 2:nNodes;
backInd = 1:nNodes-1;

flux = zeros(size(Uvec));

% rho-flux
flux(backInd) = Uvec(foreInd+varInd(2));
flux(foreInd) = flux(foreInd) - Uvec(backInd+varInd(2));

% rho*u - flux
flux(backInd+varInd(2)) = 0.5.*(3-gam).*(Uvec(foreInd+varInd(2)).^2)./Uvec(foreInd) + (gam-1).*Uvec(foreInd+varInd(3));
flux(foreInd+varInd(2)) = flux(foreInd+varInd(2))-(0.5.*(3-gam).*(Uvec(backInd+varInd(2)).^2)./Uvec(backInd) + (gam-1).*Uvec(backInd+varInd(3)));
% rho*E - flux
flux(backInd + varInd(3)) = -0.5.*(gam-1).*(Uvec(foreInd+varInd(2)).^3)./(Uvec(foreInd).^2)+gam.*Uvec(foreInd+varInd(2)).*Uvec(foreInd+varInd(3))./Uvec(foreInd);
flux(foreInd + varInd(3)) = flux(foreInd+varInd(3))-(-0.5.*(gam-1).*(Uvec(backInd+varInd(2)).^3)./(Uvec(backInd).^2)+gam.*Uvec(backInd+varInd(2)).*Uvec(backInd+varInd(3))./Uvec(backInd));

flux = flux.*0.5./dx;

end