function y = quarticTimeStep(tVec, func, y0, iterMax, matrixFlag)
    n = length(tVec);
    m = length(y0);
    y = zeros(m, n);
    dt = tVec(2) - tVec(1);
    y(:,1) = y0;
    
    A = [12*dt^2, 6*dt; 4*dt^3, 3*dt^2];
    
    for i = 1:n-1
        % Generate the quadratic spline
        yOld = [   y(1,i) + func(y(:,i), tVec(i)).*dt;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2];

        % Convert to Quartic Spline
%         coeffs = A\[yOld(1) - func(y(:,i), tVec(1))*dt - y(1,i); func(yOld, tVec(i+1)) - func(y(:,i), tVec(i))];
        if matrixFlag
            bVec = [func(yOld, tVec(i+1)) - func(y(:,i), tVec(i)); yOld(1) - func(y(:,i), tVec(1))*dt - y(1,i)];
            coeffs = A(2:-1:1, :)\bVec(end:-1:1);
        else
            a = (func(yOld, tVec(i+1)) - func(y(:,i), tVec(i)))./(4*dt^2);
            coeffs = [  a;...
                        -4*a*dt/3];
        end
        
        % Calculate next-tVec step
        yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];
        iterCt = 1;
        while norm(yNew - yOld) > 1e-5 && (iterCt <= iterMax)
            yOld = yNew;

            % Convert to Quartic Spline -> don't use matrix operations to
            % calculate this!

            if matrixFlag
                coeffs = A(2:-1:1, :)\bVec(end:-1:1);
            else
                aIter = (func(yOld, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yOld(1) - y(1,i)))./(4*dt^2);
%             coeffs = A\[yOld(1) - func(y(:,i), tVec(1))*dt - y(1,i); func(yOld, tVec(i+1)) - func(y(:,i), tVec(i))];
                coeffs = [  aIter;...
                            (-4*aIter*dt^3 + yOld(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 
            end

            % Calculate next-tVec step
            yNew = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
                        y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];
            iterCt = iterCt + 1;
        end

        y(:,i+1) = yNew;
    end

end