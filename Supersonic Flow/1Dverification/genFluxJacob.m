function fluxMat = genFluxJacob(Uvec, gam, dx)

% Inputs:
% Uvec - vector of all fluid parameters
% gam - heat coefficient ratio

% Procedure
% 1) Initialize and store nnz's into a column

nNodes = length(Uvec)/3;
varInd = (1:nNodes:length(Uvec))-1; % starting index for each eqn
foreInd = 2:nNodes;
backInd = 1:nNodes-1;

% df_1/d(rho*u)
Bcols = ones(14*(nNodes-1),1);
iVec = [backInd'; foreInd';...
        repmat([backInd'; foreInd']+varInd(2),3,1);...
        repmat([backInd'; foreInd']+varInd(3),3,1)]; % row
jVec = [[foreInd'; backInd']+varInd(2);...
        repmat([[foreInd'; backInd']; [foreInd'; backInd']+varInd(2); [foreInd'; backInd']+varInd(3)], 2,1)];% col

% df_2/drho
Bcols((2*(nNodes-1)+1):4*(nNodes-1)) = [-0.5.*(3-gam)*(Uvec((foreInd)+varInd(2))./Uvec(foreInd)).^2;...
            0.5.*(3-gam)*(Uvec((backInd)+varInd(2))./Uvec(backInd)).^2];

% df_2/d(rho*u)
Bcols((4*(nNodes-1)+1):6*(nNodes-1)) = [(3-gam).*Uvec((foreInd)+varInd(2))./Uvec(foreInd);...
            -(3-gam).*Uvec((backInd)+varInd(2))./Uvec(backInd)];

% df_2/d(rho*E)
Bcols((6*(nNodes-1)+1):7*(nNodes-1)) = gam-1;
Bcols((7*(nNodes-1)+1):8*(nNodes-1)) = 1-gam;

% df_3/d(rho)
Bcols((8*(nNodes-1)+1):10*(nNodes-1)) = [(gam-1).*(Uvec((foreInd)+varInd(2))./Uvec(foreInd)).^3 - gam.*Uvec((foreInd)+varInd(2)).*Uvec((foreInd)+varInd(3))./(Uvec(foreInd).^2);...
        -(gam-1).*(Uvec((backInd)+varInd(2))./Uvec(backInd)).^3 + gam.*Uvec((backInd)+varInd(2)).*Uvec((backInd)+varInd(3))./(Uvec(backInd).^2)];

% df_3/d(rho*u)
Bcols((10*(nNodes-1)+1):12*(nNodes-1)) = [-1.5.*(gam-1).*(Uvec((foreInd)+varInd(2))./Uvec(foreInd)).^2 + gam.*Uvec((foreInd)+varInd(3))./Uvec(foreInd);...
        1.5.*(gam-1).*(Uvec((backInd)+varInd(2))./Uvec(backInd)).^2 - gam.*Uvec((backInd)+varInd(3))./Uvec(backInd)];

% df_3/d(rho*E)
Bcols((12*(nNodes-1)+1):14*(nNodes-1)) = [gam.*Uvec((foreInd)+varInd(2))./Uvec(foreInd);...
        -gam.*Uvec((backInd)+varInd(2))./Uvec(backInd)];

% 2) Build Matrix?
fluxMat = sparse(iVec, jVec, Bcols./(2*dx));

end