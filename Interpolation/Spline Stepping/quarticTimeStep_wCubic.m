function y = quarticTimeStep_wCubic(tVec, func, y0, iterFlag)
    n = length(tVec);
    m = length(y0);
    y = zeros(m, n);
    dt = tVec(2) - tVec(1);
    y(:,1) = y0;
    
    for i = 1:n-1
        % Generate the quadratic spline
        yOld = [   y(1,i)+ func(y(:,i), tVec(i)).*dt;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2];

        % Convert to Cubic Spline
        d3y_i1 = (func(yOld, tVec(i+1)) - func(y(:,i), tVec(i)))./(6*dt);
        yCub = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*d3y_i1.*dt^2;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + d3y_i1.*dt.^3];
                
%         % Convert to Quartic Spline
%         aIter = (func(yCub, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yCub(1) - y(1,i)))./(4*dt^2);
%         coeffs = [  aIter;...
%                     (-4*aIter*dt^3 + yCub(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 
%         
%         % Calculate next-tVec step
%         yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
%                     y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];

        AA = [dt^2, dt; 12*dt^2, 6*dt];
        bVec = [(yCub(2) - 0.5.*func(y(:,i), tVec(i)).*dt^2 - y(1,i).*dt - y(2,i))./dt^2;...
                func(yCub, tVec(i+1)) - func(y(:,i), tVec(i))];
            
        coeffs1 = AA\bVec;
        
        yQuart = [  y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*coeffs1(2).*dt^2 + 4.*coeffs1(1).*dt^3;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs1(2).*dt.^3 + coeffs1(1).*dt^4];
                
        % Convert to Quartic Spline
        aIter = (func(yQuart, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yQuart(1) - y(1,i)))./(4*dt^2);
        coeffs = [  aIter;...
                    (-4*aIter*dt^3 + yQuart(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 
        
        % Calculate next-tVec step
        yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];

                
       iterCt = 1;
        while norm(yNew - yOld) > 1e-5 && iterFlag
            yOld = yNew;

%             aIter = (func(yOld, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yOld(1) - y(1,i)))./(4*dt^2);
%             coeffs = [  aIter;...
%                         (-4*aIter*dt^3 + yOld(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 

%             % Convert to Cubic Spline
%             d3y_i1 = (func(yOld, tVec(i)) - func(y(:,i), tVec(i)))./(dt);
%             yCub = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 0.5.*d3y_i1.*dt^2;...
%                         y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + (1/6)*d3y_i1.*dt.^3];
% 
%             % Convert to Quartic Spline
%             aIter = (func(yCub, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yCub(1) - y(1,i)))./(4*dt^2);
%             coeffs = [  aIter;...
%                         (-4*aIter*dt^3 + yCub(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 

%             % Calculate next-tVec step
%             yNew = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
%                         y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];

            bVec = [(yCub(2) - 0.5.*func(y(:,i), tVec(i)).*dt^2 - y(1,i).*dt - y(2,i))./dt^2;...
                    func(yOld, tVec(i+1)) - func(y(:,i), tVec(i))];

            coeffs1 = AA\bVec;

            yQuart = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*coeffs1(2).*dt^2 + 4.*coeffs1(1).*dt^3;...
                        y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs1(2).*dt.^3 + coeffs1(1).*dt^4];

            % Convert to Quartic Spline
            aIter = (func(yQuart, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yQuart(1) - y(1,i)))./(4*dt^2);
            coeffs = [  aIter;...
                        (-4*aIter*dt^3 + yQuart(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 

            % Calculate next-tVec step
            yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
                        y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];
            iterCt = iterCt + 1;
        end

        y(:,i+1) = yNew;
    end

end