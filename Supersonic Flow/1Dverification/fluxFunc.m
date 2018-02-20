function flux = fluxFunc(U, gam)

flux = [ U(2,:);...
        0.5.*(3-gam).*(U(2,:).^2)./U(1,:) + (gam-1).*U(3,:);...
        -0.5.*(gam-1).*(U(2,:).^3)./(U(1,:).^2) + gam.*U(2,:).*U(3,:)./U(1,:)];% .* g_x;

end