function threeLevel_1D(fluxFunc,visc_beta)  % three level explicit time-stepping scheme

while length(time) < 3 || norm(res(end,:)) > tol

   % Step-forward Field Values
    UU(:,:,1:2) = UU(:,:,2:3);

    % Check CFL
    if Umax*dt/dx ~= 1 %Umax*dt/dx > 1
%             error('CFL condition not met!\n');
        dt = dx./Umax;
        visc_beta = 0.5.*(VRvisc(UU(:,2:end,end),dx) + VRvisc(UU(:,1:end-1,end),dx));
        beta = max(abs(visc_beta(:)))*(dt^2)/(ratio*dx^2);
        if abs(log10(dt/dtLast)) > 0.301
            dtLast = dt;
            fprintf('Time-step changing!\nNew time step: %0.5e\nNew beta:%0.5e\n\n', dt, beta);
        end
    end

    % Calculate Next Time-Step
    UU(:,2:end-1,3) = (2.*UU(:,2:end-1,2)-(1-dt*0.5./beta)*UU(:,2:end-1,1)-(dt.^2).*flux./beta)./(1+0.5*dt./beta);
    time(end+1) = time(end)+dt;

    % Update Outflow Boundary Condition
    % 1) Extrapolate rho and E
    UU(:,end,2:3) = 5/2.*UU(:,end-1,2:3) - 2.*UU(:,end-2,2:3) + 0.5.*UU(:,end-3,2:3); % - 1/3.*UU(:,end-2,2:3);

    % 2) Fix Rho
    UU(1,end,2:3) = g_x(end).*((gam.*p_i).^(1/gam));
    
    % Update Flux and Pressure
    [flux, Umax] = fx_2Diff(@fluxFuncIsentropic, @VRvisc, UU(:,:,end), g_x, gam, dx);
    res(end+1,:) = reshape(max(abs(flux), [], 2), 1, 2);
    UU_x(:,end+1) = UU(:,xx==0.65,3);

    % Plot Residuals
    resRho.YData(end+1) = res(end,1);
    resU.YData(end+1) = res(end,2);
    Rho_x.YData(end+1) = UU_x(1,end);
    U_x.YData(end+1) = UU_x(2,end)./UU_x(1,end);
    drawnow;

end

end