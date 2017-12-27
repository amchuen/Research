function y = cubicTimeStep2(tVec, pFunc, qFunc, rFunc, y0)
    n = length(tVec);
    [m, p] = size(y0);
    y = zeros(m, p, n);
    dt = tVec(2) - tVec(1);
    y(:,:,1) = y0; % [y1' y2' ; y1 y2]
    resTol = 1e-13;
    corrTol = 1e-10;
    
    ypp = @(yMat,time) rFunc(yMat(2,:), time) - pFunc(yMat(2,:), time).*yMat(1,:) - qFunc(yMat(2,:),time).*yMat(2,:); 
%     ypp = @(yMat,time) rFunc(yMat, time) - pFunc(yMat, time).*yMat(1,:) - qFunc(yMat,time).*yMat(2,:); 

    iterCt = [0];        
    for i = 1:n-1
        
        % Build Y function at each time step -> column vector!

        % yVec = [x; y]
        func = @(yVec)  yVec.*(6 + dt^2.*qFunc(yVec,tVec(i+1)).*(1 - dt.*pFunc(yVec,tVec(i+1))./(2+dt.*pFunc(yVec,tVec(i+1)))))'...
                        - 6.*y(2,:,i)' - y(1,:,i)'.*(6*dt - dt.^2.*pFunc(yVec, tVec(i+1))./(1+0.5.*dt.*pFunc(yVec, tVec(i+1))))'...
                        - ypp(y(:,:,i),tVec(i))'.*(2*dt^2 - dt^3.*pFunc(yVec, tVec(i+1))./(2+dt*pFunc(yVec, tVec(i+1))))'...
                        - rFunc(yVec,tVec(i+1))'.*(dt^2 - dt^3.*pFunc(yVec, tVec(i+1))./(2+dt.*pFunc(yVec, tVec(i+1))))';
            
        % Perform Secant method to find root?
        if i == 1
            [yNew, iterCt(end+1)] = secantMethod(func, y(2,:,i)' - 0.25, y(2,:,i)', resTol, corrTol, 1e5);
        else
            [yNew, iterCt(end+1)] = secantMethod(func, y(2,:,i-1)', y(2,:,i)', resTol, corrTol, 1e5);
        end
        
        % Calculate first derivative
        vNew = (y(1,:,i) + 0.5.*dt.*ypp(y(:,:,i),tVec(i)) - 0.5.*dt.*qFunc(yNew,tVec(i+1)).*yNew' + 0.5.*dt*rFunc(yNew, tVec(i+1)))./(1+0.5.*dt*pFunc(yNew, tVec(i+1)));
        
        % Build new time step
        y(:,:,i+1) = [vNew; yNew'];
    end
    
    fprintf('Max no. of iterations per step: %i\n', max(iterCt));

end

function [x2, iterCt] = secantMethod(fFunc, x0, x1, fTol, xTol, maxIter)

%[xNew, res, corr, varargout]
% % Test initial guess
% xNew = xInit - jFunc(xInit)\fFunc(xInit);
% 
% % Calculate Residual and Correction
% res = norm(fFunc(xNew));
% corr = norm(xNew - xInit);

% Keep Iterating if not working

% MUST RUN INITIAL GUESS WITH LINE-SEARCH!
res = [norm(fFunc(x0)), norm(fFunc(x1))];
corr = [];
xHist = [x0, x1];

while ((res(end) > fTol) || (corr(end) > xTol)) && (length(res) <= maxIter)
    
    % Build Jacobian
    jacobian = jacobFD(x0, x1, fFunc);
    
    % Test new guess
    x2 = x1 - jacobian\fFunc(x1);

    % Calculate Residual and Correction
    res(end+1) = norm(fFunc(x2));
    corr(end+1) = norm(x2 - x1);
    xHist(:,end+1) = x2;
    
    % Re-assign
    x0 = x1;
    x1 = x2;
    
end

if (length(res) > maxIter) && (res(end) > fTol)
   fprintf('Not enough iterations to fully converge results!\n');    
end

iterCt = length(res);

end

function Jacob = jacobFD(x0, x1, func)
    Jacob = zeros(length(func(x0)), length(x0));
    for i = 1:length(x0)
        tempGuess = x0;
        tempGuess(i) = x1(i);
        tempF = func(tempGuess);
        Jacob(:,i) = (tempF - func(x0))./(x1(i) - x0(i));
    end

end