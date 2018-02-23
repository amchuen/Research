function d2Umat = matDiff2Op(dx, nNodes, nEqns)

% Notes:
% Makes 1D and 2D Difference Equations for Laplacian Operator in Matrix
% Format

% Inputs
% dx - array of different grid spacings, for each direction, 1 -> x, 2 -> y
% nNodes - no. of pts in each direction
% nEqns - number of equations

% Procedure:
% 1) make diagonal elements for 1 equation over given grid
Bvals = zeros(prod(nNodes), length(dx)+1);
Bvals(:,1) = -2.*sum((1./dx).^2);
for i = 1:length(dx)
    Bvals(:,i+1) = (1./dx(i)).^2;
    if i == 1
        Bvals(nNodes(i)+1:nNodes(i):prod(nNodes),i+1) = 0;
    end
end

% 2) extend for number of equations
Bvals = repmat(Bvals, nEqns, 1);

% 3) build sparse matrix
dimArr = [0, 1];
if length(dx) == 2
    dimArr(end+1) = nNodes(1);
end    

d2Umat = spdiags(Bvals, dimArr, size(Bvals,1), size(Bvals,1));
d2Umat = d2Umat + d2Umat';

end