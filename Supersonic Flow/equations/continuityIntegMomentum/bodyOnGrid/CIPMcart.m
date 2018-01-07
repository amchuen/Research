function [fvOUT, waveSpd] = CIPMcart(GR, FL, BC, FF)

fvOUT = zeros(size(FF));

%% Calculate Operations on Continuity

% Density averages
rho_e = 0.5.*(FF(2:end-1,3:end,1) + FF(2:end-1,2:end-1,1));
rho_w = 0.5.*(FF(2:end-1,1:end-2,1) + FF(2:end-1,2:end-1,1));
rho_n = 0.5.*(FF(2:end-1,2:end-1,1) + FF(3:end,2:end-1,1));
rho_s = 0.5.*(FF(2:end-1,2:end-1,1) + FF(1:end-2,2:end-1,1));

% Velocity -> first order
phiX_w = diff(FF(2:end-1,1:end-1,2),1,2)./GR.dx;
phiX_e = diff(FF(2:end-1,2:end,2),1,2)./GR.dx;
phiY_n = diff(FF(2:end,2:end-1,2),1,1)./GR.dy;
phiY_s = diff(FF(1:end-1,2:end-1,2),1,1)./GR.dy;

fvOUT(2:end-1,2:end-1,1) = (rho_e.*phiX_e - rho_w.*phiX_w)./GR.dx + (rho_n.*phiY_n - rho_s.*phiY_s)./GR.dy;

%% Calculate Operations on Potential
phiX = 0.5.*(phiX_w + phiX_e);
phiY = 0.5.*(phiY_n + phiY_s);
PP = (FF(2:end-1,2:end-1,1).^(FL.gam-1) - 1)./((FL.gam - 1).*FL.M0^2);
fvOUT(2:end-1,2:end-1,2) = PP + 0.5.*((phiX).^2 + (phiY).^2 - 1);

%% Calculate Stability Output
rho = FF(:,:,1);
waveSpd = [max([max(abs(phiX(:))), 1, max(rho(:))]), max([max(abs(phiY(:))), 1, max(rho(:))])];

end