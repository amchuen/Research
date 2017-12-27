function y = cubicTimeStep(tVec, func, y0, iterFlag)
    n = length(tVec);
    [m, p] = size(y0);
    y = zeros(m, p, n);
    dt = tVec(2) - tVec(1);
    y(:,:,1) = y0;
    
    %             yNew = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 0.5.*d3y_i1.*dt^2;...
    %                         y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + (1/6)*d3y_i1.*dt.^3];
    Amat = [1 0 dt 3*dt^2; dt 1 0.5.*dt^2 dt^3];
    
    for i = 1:n-1
        d_j = zeros(1,p);
        % Generate the quadratic spline
        yOld = Amat * [y(:,:,i); func(y(:,:,i), tVec(i)); d_j];

        % Convert to Cubic Spline
        d_j = (func(yOld, tVec(i+1)) - func(y(:,:,i), tVec(i)))./(6*dt);

%         Calculate next-tVec step
        yNew = Amat * [y(:,:,i); func(y(:,:,i), tVec(i)); d_j];
%         while norm(yNew - yOld) > 1e-5 && iterFlag
        if iterFlag
            for ii = 1:3
                yOld = yNew;

                % Convert to Cubic Spline
                d_j = (func(yOld, tVec(i+1)) - func(y(:,:,i), tVec(i)))./(6*dt);

                % Calculate next-tVec step
                yNew = Amat * [y(:,:,i); func(y(:,:,i), tVec(i)); d_j];
            end
        end

        y(:,:,i+1) = yNew(1:2,:);
    end

end