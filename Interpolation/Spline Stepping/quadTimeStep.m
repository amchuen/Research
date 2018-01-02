function y = quadTimeStep(tVec, func, y0, iterFlag)
    n = length(tVec);
    [m, p] = size(y0);
    y = zeros(m, p, n);
    dt = tVec(2) - tVec(1);
    y(:,:,1) = y0;
    
    for i = 1:n-1
        % Generate the quadratic spline
        yNew = [    y(1,:,i) + func(y(:,:,i), tVec(i)).*dt;...
                    y(2,:,i) + y(1,:,i).*dt + 0.5.*func(y(:,:,i), tVec(i)).*dt^2];

        % Convert to Cubic Spline
%         d3y_i1 = (func(yOld, tVec(i+1)) - func(y(:,:,i), tVec(i)))./(dt);
        y(:,:,i+1) = yNew;
    end

end