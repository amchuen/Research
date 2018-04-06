function flux = fluxFuncIsentropic(U, gam)

flux = [ U(2,:);...
        (U(2,:).^2)./U(1,:) + (U(1,:).^gam)./gam];% .* g_x;

end