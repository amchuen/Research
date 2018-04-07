function [FF, PP] = fluxFuncIsentropic(U, gam)

% Nonlinear Transport Terms
FF = [ U(2,:);...
        (U(2,:).^2)./U(1,:) + (U(1,:).^gam)./gam];% .* g_x;
    
% Pressure
PP = (U(1,:).^gam)./gam;

end