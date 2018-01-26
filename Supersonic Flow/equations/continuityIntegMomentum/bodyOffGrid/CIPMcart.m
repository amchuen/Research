function fvOUT = CIPMcart(GR, FL, BC, FF)

fvOUT = zeros(size(FF));

%% Calculate Operations on Continuity

% Density averages
rho_e = 0.5.*[(FF(:,2:end,1) + FF(:,1:(end-1),1)), (FF(:,end,1) + bcCalc(GR, BC, FF, 'E', 1))];
rho_w = 0.5.*[(bcCalc(GR, BC, FF, 'W', 1) + FF(:,1,1)), (FF(:,1:(end-1),1) + FF(:,2:end,1))];
rho_n = 0.5.*[(FF(2:end,:,1) + FF(1:(end-1),:,1)); (bcCalc(GR, BC, FF, 'N', 1) + FF(end,:,1))];
rho_s = 0.5.*[bcCalc(GR, BC, FF, 'S', 1)+FF(1,:,1); FF(2:end,:,1) + FF(1:(end-1),:,1)];

% Velocity -> first order
phiX_w = diff([bcCalc(GR, BC, FF, 'W', 2), FF(:,:,2)],1,2)./GR.dx;
phiX_e = diff([FF(:,:,2), bcCalc(GR, BC, FF, 'E', 2)],1,2)./GR.dx;
phiY_n = diff([FF(:,:,2); bcCalc(GR, BC, FF, 'N', 2)],1,1)./GR.dy;
phiY_s = diff([bcCalc(GR, BC, FF, 'S', 2); FF(:,:,2)],1,1)./GR.dy;

fvOUT(:,:,1) = (rho_e.*phiX_e - rho_w.*phiX_w)./GR.dx + (rho_n.*phiY_n - rho_s.*phiY_s)./GR.dy;

%% Calculate Operations on Potential

PP = (FF(:,:,1).^(FL.gam) - 1)./((FL.gam - 1).*FL.M0^2);
fvOUT(:,:,2) = PP + 0.5.*((0.5.*(phiX_w + phiX_e)).^2 + (0.5.*(phiY_n + phiY_s)).^2 - 1);

end