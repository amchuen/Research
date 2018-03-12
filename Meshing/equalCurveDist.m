function xLoc = equalCurveDist(coeffs, xRange, npts)

dS = @(x) sqrt(1 + (polyval(((length(coeffs)-1):-1:1)'.*coeffs(1:end-1), x)).^2);
totLen = integral(dS, xRange(1), xRange(2));

lenDist = linspace(0, totLen, npts);
xLoc = zeros(size(lenDist)); xLoc(1) = xRange(1); xLoc(end) = xRange(2);
for i = 2:(npts-1)
    lenFunc = @(x) integral(dS, xRange(1), x) - lenDist(i);
    xLoc(i) = newtonSys(lenFunc, dS, xLoc(i), 1e-6, 1e-6, 1e5, 0); 
end

end

function varargout = newtonSys(fFunc, jFunc, xNew, fTol, xTol, maxIter, lsFlag)

%[xNew, res, corr, varargout]
% % Test initial guess
% xNew = xInit - jFunc(xInit)\fFunc(xInit);
% 
% % Calculate Residual and Correction
% res = norm(fFunc(xNew));
% corr = norm(xNew - xInit);

% Keep Iterating if not working

% MUST RUN INITIAL GUESS WITH LINE-SEARCH!
xOld = xNew;
xNew = xOld - jFunc(xOld)\fFunc(xOld);
res = norm(fFunc(xNew));
corr = norm(xNew - xOld);
% xHist = xNew;

while ((res(end) > fTol) || (corr(end) > xTol)) && (length(res) <= maxIter)
    % Test new guess
    xOld = xNew;
    xNew = xOld - jFunc(xOld)\fFunc(xOld);

    % Line Search Option
    if lsFlag
        lambda = 1;
        xLS = xOld + lambda*(xNew-xOld);
        nlambdas = 1;
%         alpha=1e-4;
        while (norm(fFunc(xLS))^2 > norm(fFunc(xOld))^2) && (nlambdas <= maxIter)
           lambda = lambda * 0.5;
           xLS = xOld + lambda*(xNew-xOld);
           nlambdas = nlambdas + 1;
        end
       
        xNew = xLS;
        
    end
    
    % Calculate Residual and Correction
    res(end+1) = norm(fFunc(xNew));
    corr(end+1) = norm(xNew - xOld);
%     xHist(:,end+1) = xNew;    
    
end

if (length(res) > maxIter) && (res(end) > fTol)
   fprintf('Not enough iterations to fully converge results!\n');    
end

varargout{1} = xNew;
varargout{2} = res;
varargout{3} = corr;
% varargout{4} = xHist;

end