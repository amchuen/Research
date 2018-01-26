function [fvOUT, waveSpd] = CIPMcyl(GR, FL, BC, FF)

fvOUT = zeros(size(FF));

%% Calculate Operations on Continuity

% Density averages
rho_e = 0.5.*[(FF(:,2:end,1) + FF(:,1:(end-1),1)), (FF(:,end,1) + bcCalc(GR, BC, FF, 'E', 1))];
rho_w = 0.5.*[(bcCalc(GR, BC, FF, 'W', 1) + FF(:,1,1)), (FF(:,1:(end-1),1) + FF(:,2:end,1))];
rho_n = 0.5.*[(FF(2:end,:,1) + FF(1:(end-1),:,1)); (bcCalc(GR, BC, FF, 'N', 1) + FF(end,:,1))];
rho_s = 0.5.*[bcCalc(GR, BC, FF, 'S', 1)+FF(1,:,1); FF(2:end,:,1) + FF(1:(end-1),:,1)];

% Velocity -> first order
if GR.isPolar
    phi2_w = diff([bcCalc(GR, BC, FF, 'W', 2), FF(:,:,2)],1,2)./(GR.RR.*GR.dT);
    phi2_e = diff([FF(:,:,2), bcCalc(GR, BC, FF, 'E', 2)],1,2)./(GR.RR.*GR.dT);
    phi1_n = diff([FF(:,:,2); bcCalc(GR, BC, FF, 'N', 2)],1,1)./GR.dR;
    phi1_s = diff([bcCalc(GR, BC, FF, 'S', 2); FF(:,:,2)],1,1)./GR.dR; 
    
    fvOUT(:,:,1) = (rho_e.*phi2_e - rho_w.*phi2_w)./(GR.dT.*GR.RR) + (GR.RR_N.*rho_n.*phi1_n - GR.RR_S.*rho_s.*phi1_s)./(GR.dR .* GR.RR);
else
    phi2_w = diff([bcCalc(GR, BC, FF, 'W', 2), FF(:,:,2)],1,2)./GR.dx;
    phi2_e = diff([FF(:,:,2), bcCalc(GR, BC, FF, 'E', 2)],1,2)./GR.dx;
    phi1_n = diff([FF(:,:,2); bcCalc(GR, BC, FF, 'N', 2)],1,1)./GR.dy;
    phi1_s = diff([bcCalc(GR, BC, FF, 'S', 2); FF(:,:,2)],1,1)./GR.dy;
    
    fvOUT(:,:,1) = (rho_e.*phi2_e - rho_w.*phi2_w)./GR.dx + (rho_n.*phi1_n - rho_s.*phi1_s)./GR.dy;
end


%% Calculate Operations on Potential
phi2 = 0.5.*(phi2_w + phi2_e);
phi1 = 0.5.*(phi1_n + phi1_s);
PP = (FF(:,:,1).^(FL.gam-1) - 1)./((FL.gam - 1).*FL.M0^2);
fvOUT(:,:,2) = PP + 0.5.*(phi2.^2 + phi1.^2 - 1);

%% Calculate Stability Output
% rho = FF(:,:,1);
phi2w = phi2;%./(GR.RR.*2.*GR.dT) ;%./GR.RR;
phi1w = phi1;%./(2.*GR.dR);%./GR.RR;
% rhow = rho;%./GR.RR;
% waveSpd = [max([max(abs(phi2w(:))), max(rhow(:))]), max([max(abs(phi1w(:))), max(rhow(:))])];
waveSpd = [max(abs(phi2w(:))), max(abs(phi1w(:)))];

end