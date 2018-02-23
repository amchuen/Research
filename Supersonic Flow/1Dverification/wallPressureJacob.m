function Pmat = wallPressureJacob(Uvec, dgdx, gam)

Bvals = zeros(size(Uvec));
nNodes = length(Uvec)/3;
iVec = repmat(1:nNodes, 1,3) + nNodes; 
jVec = 1:length(Uvec);

% Initialize Values
Bvals(1:nNodes) = dgdx.*(gam-1).*(Uvec((1:nNodes)+2*nNodes) + 0.5.*(Uvec((1:nNodes)+nNodes)./Uvec(1:nNodes)).^2);
Bvals((1:nNodes)+nNodes) = dgdx.*(gam-1).*(Uvec((1:nNodes)+2*nNodes) - Uvec((1:nNodes)+nNodes)./Uvec(1:nNodes));
Bvals((1:nNodes)+2*nNodes) = dgdx.*(gam-1);

% Build Matrix
Pmat = sparse(iVec, jVec, Bvals, length(Uvec), length(Uvec));

end