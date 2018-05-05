function [flux, PP, varargout] = flux_FullEuler(U, gam)

%% Calculate Nonlinear Terms
flux = [ U(2,:);...
        0.5.*(3-gam).*(U(2,:).^2)./U(1,:) + (gam-1).*U(3,:);...
        -0.5.*(gam-1).*(U(2,:).^3)./(U(1,:).^2) + gam.*U(2,:).*U(3,:)./U(1,:)];% .* g_x;
    
%% Pressure Calculation
PP = (gam-1).*(U(3,:) - 0.5.*(U(2,:).^2) ./ U(1,:));

%% Characteristic Speeds
CC = sqrt(gam.*PP./U(1,:));
vel = U(2,:)./U(1,:);

%% Misc. Outputs

% characteristic speeds
varargout{1} = cat(1, vel, vel+CC, vel-CC);

end