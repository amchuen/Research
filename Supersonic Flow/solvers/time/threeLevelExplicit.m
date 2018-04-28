function OUT = threeLevelExplicit(GR, FL, BC, U0, fluxFunc, varargin)

%% Run Checks

% Generic three-level time-stepper 

% Check if it is fixed time-sim or run until steady state

% Check if cylindrical or cartesian coordinates
% if GR.isPolar
%     % Raidus Averagess
%     GR.RR_N = [0.5.*(GR.RR(2:end,:) + GR.RR(1:(end-1),:)); 0.5.*((GR.RR(end,:)+GR.dr) + GR.RR(end,:))];
%     GR.RR_S = 0.5.*([2.*GR.RR(1,:)-GR.dr; GR.RR(2:end,:) + GR.RR(1:(end-1),:)]);
% end

%% Initialize Variables

% Matrix Dimensions
% 1:Y, 2:X, 3:vec, 4:t
UU = repmat(U0,1,1,1,3);%struct('fv',repmat(U0,1,1,1,3),'f2',zeros(size(U0)), 'f1', zeros(size(U0)));

% Calculate Residual before running simulation
res = ones(2,size(UU,3)); %resCalc(GR, FL, BC, func, epsFunc, UU(:,:,:,end));
maxRes = max(res,[],1);
time = [-GR.dt, 0];

%% Initialize Plots for Residual
figure(1);
resPlot = semilogy(res);
legend(BC.N.varName, 'Location', 'Best');
title(['t = ' num2str(time(end))]);

%% Run Simulation
dtLast = GR.dt;
while norm(res(end,:)./maxRes) > GR.tol || time(end) < GR.tEnd
    
    % Step-Forward
    UU(:,:,:,1:2) = UU(:,:,:,2:3);
    
    % Calculate Flux Terms
    [flux, waveSpd, ~, eps1, eps2] = fluxFunc(GR, FL, BC, UU(:,:,:,2));
        
    % Check CFL
    cflWave = GR.dt.*(waveSpd(1)./GR.dx + waveSpd(2)./GR.dy);
    GR.cflFactor = 1;
    if (max(abs(cflWave(:))) ~= GR.CFL)% -> DF method might allow for unconditional stability for pure diffusion
        GR.cflFactor = GR.CFL./max(cflWave(:));
        GR.dt = GR.dt .* GR.cflFactor;
        
        if abs(log10(GR.dt/dtLast)) >= 0.5
            fprintf('CFL condition not met!\n');
            fprintf('Changing time steps!\n');
            fprintf('New time step:%0.5e\n', GR.dt);
            dtLast = GR.dt;
        end
    end
    
    % Compute Time-Step
    beta = (GR.dt^2).*(eps1./(GR.dx^2) + eps2./(GR.dy^2))./GR.ratio;
    if GR.isPolar
        UU(:,:,:,3) = ((1-alpha2-alpha1).*UU(:,:,:,1)+alpha2.*(UU.f2)+alpha1.*(UU.f1)./(0.5.*(GR.RR_N+GR.RR_S))-(1+1/cflFactor).*GR.dt.*(flux - epsFunc(GR,BC,'X').*rotLaplace))./(1+alpha2+alpha1);
    else
        UU(:,:,:,3) = (2.*beta.*UU(:,:,:,2)./(GR.dt^2)-(beta./(GR.dt^2)-0.5./GR.dt).*UU(:,:,:,1) - flux)./(beta./(GR.dt^2)+0.5./GR.dt);
    end
    
    % Calculate Residual
    res(end+1,:) = resCalc(GR, UU(:,:,:,end-1:end));    
    time(end+1) = time(end) + GR.dt;
    
    runChecks(UU(:,:,:,end));
    
    if (length(res) >= 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', norm(res(end,:)));
    end
    
    % Plot Residuals
    for iii = 1:size(res,2)
        resPlot(iii).YData(end+1) = res(end,iii);
    end
    drawnow;
    
    if size(res,1) < 5
        maxRes = max(res, [], 1);
    else
        maxRes = max(res(1:5,:), [], 1);
    end
    
end


%% Output for Post-Processing
OUT.Uvals = UU;
OUT.time = time;
OUT.res = res;

end

function runChecks(FV)

    % Run Additional Time-Stepping Checks
    if any(any(FV(:,:,1,end)<=0))
        fprintf('Check Rho!\n');
    end
    
    if ~isreal(FV(:))% || any(isnan(FV(:))))
        error('Solution exhibits non-real values.\n');
%         break;
    end
    
    if  any(isnan(FV(:)))
        error('Solution exhibits NaN at nodes. Check if divided by zero!.\n');

    end
    
    if ~isreal(FV(:))
       fprintf('Check phi!\n'); 
    end
    
end

function res = resCalc(GR, Uvals) 

% checks steady state by calculateing Ut?
for i = 1:size(Uvals,3)
    tempres = abs(Uvals(:,:,i,end) - Uvals(:,:,i,end-1))./GR.dt;
    res(i) = max(tempres(:));
end

end