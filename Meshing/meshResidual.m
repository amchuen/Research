function res = meshResidual(UU)

% Get Derivatives?
xZeta = UU(2:end-1, 3:end, 1) - UU(2:end-1, 1:end-2,1);
yZeta = UU(2:end-1, 3:end, 2) - UU(2:end-1, 1:end-2,2);
xEta = UU(3:end,2:end-1,1) - UU(1:end-2,2:end-1,1);
yEta = UU(3:end,2:end-1,2) - UU(1:end-2,2:end-1,2);

% Calculate alpha, beta, gamma, jacobian
alpha = xEta.^2 + yEta.^2;
beta = xZeta.*xEta + yZeta.*yEta;
gamma = xZeta.^2 + yZeta.^2;
jac = xZeta.*yEta - yZeta.*xEta;

% Calculate 2nd Derivs
U_ZZ = UU(2:end-1, 3:end,:) - 2.*UU(2:end-1, 2:end-1,:) + UU(2:end-1, 1:end-2,:);
U_EE = UU(3:end, 2:end-1,:) - 2.*UU(2:end-1, 2:end-1,:) + UU(1:end-2, 2:end-1,:);
U_ZE = UU(3:end, 3:end,:) - UU(1:end-2, 3:end,:) - UU(3:end, 1:end-2, :) + UU(1:end-2, 1:end-2, :);
FF = cat(3, (jac.^2).*(UU(2:end-1, 2:end-1,3).*xZeta + UU(2:end-1, 2:end-1, 4).*xEta),...
            (jac.^2).*(UU(2:end-1, 2:end-1,3).*yZeta + UU(2:end-1, 2:end-1, 4).*yEta),...
            zeros(size(UU,1)-2, size(UU,2)-2, 2));
errU = alpha.*U_ZZ + gamma.*U_EE -2.*beta.*U_ZE + FF;

res = max(abs(errU(:)));

end