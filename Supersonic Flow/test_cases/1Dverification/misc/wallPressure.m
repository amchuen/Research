function Pvec = wallPressure(Uvec, dgdx, gam)

Pvec = zeros(size(Uvec));
nNodes = length(Uvec)/3;

Pvec((1:nNodes)+nNodes) = dgdx.*(gam-1).*(Uvec((1:nNodes)+2*nNodes) - 0.5.*(Uvec((1:nNodes)+nNodes).^2)./Uvec(1:nNodes));

end