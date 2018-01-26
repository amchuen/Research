function eulerIsentropic (GR, FL, BC, CT, RHO, AA, BB, PP)

res = [1, 1, 1];
tic;
while (max(res(end, :)) > tol*max(res(:)))|| (size(res,1) < CT.iter_min) % iterate through time
    
    % Update Field Values
    RHO.fv(:,:,1:2) = RHO.fv(:,:,2:3);
    AA.fv(:,:,1:2) = AA.fv(:,:,2:3);
    BB.fv(:,:,1:2) = BB.fv(:,:,2:3);
    PP.fv = RHO.fv(:,:,2).^(FL.gam) ./ (FL.gam .* FL.M0.^2); % pressure
    
    % Update Boundary Conditions
%     BCval_BB = {0, 0, 0, RHO.fv(1,:,2).*dyBdx};
%     BCval_AB_R = {0, 0, 0, AA.fv(1,:, 2).*dyBdx};
%     BB.BC = init_fv(BC_BB, BCval_BB, GR.XX);
%     AB_R.BC = init_fv(BC_AB_R, BCval_AB_R, GR.XX);
    BB.BC.S = RHO.fv(1,:,2).*dyBdx;
    AB_R.BC.S = AA.fv(1,:, 2).*dyBdx;
    
    % Check CFL conditions
    Ux = (AA.fv(:,:,2)./RHO.fv(:,:,2));
    Vy = (BB.fv(:,:,2)./RHO.fv(:,:,2));
    CFL_i = (max(abs(Ux(:)))./dx + max(abs(Vy(:)))./dy)*dt;

    if CFL_i >= 1.0 && CFL_on
       fprintf('CFL condition not met!\n');
%            if CFL_on
       fprintf('Decreasing time steps!\n');
       dt = dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
%            end
%     elseif CFL_i < 0.8
%         dt = dt / CFL_i;
    end
    
    % Calculate new RHO
    [Ax, ~, ~] = grad_f(AA.fv(:,:,2), 2, dx, AA.BC, 1);
    [By, ~, ~] = grad_f(BB.fv(:,:,2), 1, dy, BB.BC, 1);
	laplace_RHO = laplace_f(RHO.fv(:,:,2), dx, dy, RHO.BC, 1);
    RHO.fv(:,:,3) = (eps_s.*laplace_RHO - Ax - By - (eps_t./(dt^2) - 0.5./dt).*RHO.fv(:,:,1) + 2.*eps_t.*RHO.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
        
    % Calculate new A
    [A2_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).^2)./RHO.fv(:,:,2), 2, dx, AA.BC, 1);
    [AB_rho_y, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 1, dy, AB_R.BC, 1);
    laplace_A = laplace_f(AA.fv(:,:,2), dx, dy, AA.BC, 1);
    [P_x, ~, ~] = grad_f(PP.fv, 2, dx, PP.BC, 1);
    AA.fv(:,:,3) = (eps_s.*laplace_A - P_x - A2_rho_x - AB_rho_y - (eps_t./(dt^2) - 0.5./dt).*AA.fv(:,:,1) + 2.*eps_t.*AA.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    % Calculate new B
    [AB_rho_x, ~, ~] = grad_f((AA.fv(:,:,2).*BB.fv(:,:,2))./RHO.fv(:,:,2), 2, dx, AB_R.BC, 1);
    [B2_rho_y, ~, ~] = grad_f((BB.fv(:,:,2).^2)./RHO.fv(:,:,2), 1, dy, BB.BC, 1);
    laplace_B = laplace_f(BB.fv(:,:,2), dx, dy, BB.BC, 1);
    [P_y, ~, ~] = grad_f(PP.fv, 1, dy, PP.BC, 1);
    BB.fv(:,:,3) = (eps_s.*laplace_B - P_y - AB_rho_x - B2_rho_y - (eps_t./(dt^2) - 0.5./dt).*BB.fv(:,:,1) + 2.*eps_t.*BB.fv(:,:,2)./(dt.^2))./(eps_t./(dt^2) + 0.5./dt);
    
    % Enforce Dirichlet BC's
    RHO.fv(:,1,3) = RHO.BC.W;
    RHO.fv(end,:,3) = RHO.BC.N;
    AA.fv(:,1,3) = AA.BC.W;
    AA.fv(end,:,3) = AA.BC.N; 
    BB.fv(:,1,3) = BB.BC.W;
    BB.fv(end,:,3) = BB.BC.N;
    BB.fv(1,:,3) = BB.BC.S;
    
    % Calculate Residuals
    R_err = (abs(RHO.fv(:,:,3) - RHO.fv(:,:,2)));
    A_err = (abs(AA.fv(:,:,3) - AA.fv(:,:,2)));
    B_err = (abs(BB.fv(:,:,3) - BB.fv(:,:,2)));
    
    if (size(res,1) == 1) && all(res(end,:) == 1)
        res(1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))]; 
    else
        res(end+1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))];
    end
    
    if (size(res, 1) > 500) && (mod(size(res, 1), 2000) == 0)
        fprintf('Iteration Ct: %i\n', size(res, 1));
        fprintf('Current Residual: %0.5e\n', max(res(end, :)));
        toc;
        figure(1);semilogy(1:size(res,1), res(:,1), 1:size(res,1), res(:,2), 1:size(res,1), res(:,3));
        hold off;
        legend('Density', '\rho u', '\rho v');
        fprintf('\n');
    end
    
    if ~isreal(RHO.fv(:,:,end)) || any(any(isnan(RHO.fv(:,:,end))))
        fprintf('Density exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    if ~isreal(AA.fv(:,:,end)) || any(any(isnan(AA.fv(:,:,end))))
        fprintf('Vx exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
    if ~isreal(BB.fv(:,:,end)) || any(any(isnan(BB.fv(:,:,end))))
        fprintf('Vy exhibits non-solutions (either non-real or NaN) in nodes!\n');
        break;
    end
    
end

%% OUTPUT

Ux = (AA.fv(:,:,2)./RHO.fv(:,:,2));
Vy = (BB.fv(:,:,2)./RHO.fv(:,:,2));

q2_ij = (Ux).^2 + (Vy).^2;

OUT.RHO = RHO;
OUT.AA = AA;
OUT.BB = BB;

    
    
end