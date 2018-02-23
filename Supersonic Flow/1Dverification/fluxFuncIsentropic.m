function flux = fluxFuncIsentropic(U, gam)

flux = [ U(2,:);...
        (U(2,:).^2)./U(1,:) + (3.*U(1,:) - 0.5.*(U(2,:).^2)./U(1,:)).*(gam-1)./gam];% .* g_x;

end