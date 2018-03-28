function OUT = dufortFrankel(GR, FL, BC, func, epsFunc, U0)

%% Run Checks

% Check if it is fixed time-sim or run until steady state

% Check if cylindrical or cartesian coordinates

%% Initialize Variables

% Matrix Dimensions
% 1:Y, 2:X, 3:vec, 4:t
UU = struct('fv',repmat(U0,1,1,1,2),'fx',zeros(size(U0)), 'fy', zeros(size(U0)));

% Calculate Residual before running simulation
res = [1, 1; 1,1]; %resCalc(GR, FL, BC, func, epsFunc, UU.fv(:,:,:,end));
time = [-GR.dt, 0];

%% Run Simulation

while norm(res(end,:)) > GR.tol || time(end) < GR.tEnd
    
    % Define stability coefficients
    alphaX = 2.*GR.dt.*epsFunc(GR, BC, 'X')./(GR.dx.^2);
    alphaY = 2.*GR.dt.*epsFunc(GR, BC, 'Y')./(GR.dy.^2);
    alphaTot = alphaX + alphaY;
    
    % Check CFL
    cflFactor = 1;
    if any(alphaTot(:) ~= GR.CFL)
        cflFactor = GR.CFL./max(alphaTot(:));
        alphaX = alphaX .* (1+cflFactor)./2;
        alphaY = alphaY .* (1+cflFactor)./2;
        GR.dt = GR.dt .* cflFactor;
        
        fprintf('CFL condition not met!\n');
        fprintf('Decreasing time steps!\n');
        fprintf('New time step:%0.5e\n', GR.dt);
    end
    
    % Compute Center-Spacing Derivatives
    UU.fx(:,2:end-1,:) = UU.fv(:,3:end,:,end) + UU.fv(:,1:end-2,:,end);
    UU.fx(:,1,:) = UU.fv(:,2,:,end) + bcCalc(GR, BC, UU.fv(:,:,:,end), 'W');
    UU.fx(:,end,:) = UU.fv(:,end-1,:,end) + bcCalc(GR, BC, UU.fv(:,:,:,end), 'E');
    UU.fy(2:end-1,:,:) = UU.fv(3:end,:,:,end) + UU.fv(1:end-2,:,:,end);
    UU.fy(1,:,:) = UU.fv(2,:,:,end) + bcCalc(GR, BC, UU.fv(:,:,:,end),'S');
    UU.fy(end,:,:) = UU.fv(end-1,:,:,end) + bcCalc(GR,BC, UU.fv(:,:,:,end),'N');
    
    % Compute Time-step
    UU.fv(:,:,:,end+1) = ((1-alphaX-alphaY).*UU.fv(:,:,:,end-1)+alphaX.*(UU.fx)+alphaY.*(UU.fy)-(1+1/cflFactor).*GR.dt.*func(GR, FL, BC, UU.fv(:,:,:,end)))./(1+alphaX+alphaY);
    
    % Calculate Residual
    res(end+1,:) = resCalc(GR, UU.fv(:,:,:,end:-1:end-2));    
    time(end+1) = time(end) + GR.dt;
    
    % Run Additional Time-Stepping Checks
    if any(any(UU.fv(:,:,1,end)<=0))
        fprintf('Check Rho!\n');
    end
    
    if (~isreal(UU.fv(:,:,:,end)) || any(any(any(isnan(UU.fv(:,:,:,end))))))
        error('Solution exhibits non-solutions (either non-real or NaN) in nodes!\n');
%         break;
    end
    
    if ~isreal(UU.fv(:,:,:,end))
       fprintf('Check phi!\n'); 
    end
    
    if (length(res) >= 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
    end
    
    figure(1);
    semilogy(time, res(:,1));
    hold on; semilogy(time, res(:,2));
    legend('\rho', '\phi', 'Location', 'Best');
    %contourf(GR.XX, GR.YY, UU.fv(:,:,1,end),50);
    title(['t = ' num2str(time(end))]);
    hold off;
    %colorbar;
    drawnow;
    
end


%% Output for Post-Processing
OUT.Uvals = UU.fv;
OUT.time = time;
OUT.res = res;

end

function res = resCalc(GR, Uvals) 

% checks steady state by calculateing Ut?
for i = 1:size(Uvals,3)
    tempres = abs(Uvals(:,:,i,end) - Uvals(:,:,i,end-1))./GR.dt;
    res(i) = max(tempres(:));
end

end