function y = quarticTimeStep_iterQuart(tVec, func, y0, iterFlag)
    n = length(tVec);
    m = length(y0);
    y = zeros(m, n);
    dt = tVec(2) - tVec(1);
    y(:,1) = y0;
    
    A1 = diag([1,1,2,0,0]); A1(4,:) = [0 0 2 6*dt 12*dt^2];
    A2 = A1;
    A1(end,:) = dt.^(0:4);
    A2(end,:) = (0:4).*dt.^(-1:3);
    
    Anew = [A2(end,:); A1(end,:)];
     
    for i = 1:n-1
        % Generate the quadratic spline
        yOld = [   y(1,i)+ func(y(:,i), tVec(i)).*dt;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2];

        % Convert to Cubic Spline
        d3y_i1 = (func(yOld, tVec(i+1)) - func(y(:,i), tVec(i)))./(6*dt);
        yCub = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*d3y_i1.*dt^2;...
                    y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + d3y_i1.*dt.^3];
                
        % Calculate new ydot using previous step, yddot, and y
        b1 = [y(end:-1:1,i); func(y(:,i), tVec(i)); func(yCub, tVec(i+1)); yCub(end)];
        coeffs1 = A1 \ b1;
        
        yQuart = Anew * coeffs1;
        
        b2 = b1;
        b2(end-1:end) = [func(yQuart, tVec(i+1)); yQuart(1)];
        coeffs2 = A2 \ b2;
        
        yNew = Anew * coeffs2;
        
        
%         AA = [dt^2, dt; 12*dt^2, 6*dt];
%         bVec = [(yCub(2) - 0.5.*func(y(:,i), tVec(i)).*dt^2 - y(1,i).*dt - y(2,i))./dt^2;...
%                 (func(yCub, tVec(i+1)) - func(y(:,i), tVec(i)))];
%             
%         coeffs1 = AA\bVec;
%         
%         yQuart = [  y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*coeffs1(2).*dt^2 + 4.*coeffs1(1).*dt^3;...
%                     y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs1(2).*dt.^3 + coeffs1(1).*dt^4];
%                 
%         % Calculate new y using previous step, ydot, and yddot
%         aIter = (func(yQuart, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yQuart(1) - y(1,i)))./(4*dt^2);
%         coeffs = [  aIter;...
%                     (-4*aIter*dt^3 + yQuart(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 
%         
%         % Calculate next-tVec step
%         yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
%                     y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];

                
       iterCt = 1;
        while norm(yNew - yOld) > 1e-5 && iterFlag
            yOld = yNew;
            
            b1 = [y(end:-1:1,i); func(y(:,i), tVec(i)); func(yOld, tVec(i+1)); yOld(end)];
            coeffs1 = A1 \ b1;

            yQuart = Anew * coeffs1;

            b2 = b1;
            b2(end-1:end) = [func(yQuart, tVec(i+1)); yQuart(1)];
            coeffs2 = A2 \ b2;

            yNew = Anew * coeffs2;

%             bVec = [(yCub(2) - 0.5.*func(y(:,i), tVec(i)).*dt^2 - y(1,i).*dt - y(2,i))./dt^2;...
%                     func(yOld, tVec(i+1)) - func(y(:,i), tVec(i))];
% 
%             coeffs1 = AA\bVec;
% 
%             yQuart = [    y(1,i)+ func(y(:,i), tVec(i)).*dt + 3.*coeffs1(2).*dt^2 + 4.*coeffs1(1).*dt^3;...
%                         y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs1(2).*dt.^3 + coeffs1(1).*dt^4];
% 
%             % Convert to Quartic Spline
%             aIter = (func(yQuart, tVec(i+1)) + func(y(:,i), tVec(i)) - 2/dt * (yQuart(1) - y(1,i)))./(4*dt^2);
%             coeffs = [  aIter;...
%                         (-4*aIter*dt^3 + yQuart(1) - y(1,i) - func(y(:,i),tVec(i)).*dt)/(3*dt^2)]; 
% 
%             % Calculate next-tVec step
%             yNew = [    y(1,i) + func(y(:,i), tVec(i)).*dt + 3.*dt^2.*coeffs(2) + 4.*coeffs(1).*dt^3;...
%                         y(2,i) + y(1,i).*dt + 0.5.*func(y(:,i), tVec(i)).*dt^2 + coeffs(2).*dt.^3 + coeffs(1).*dt^4];
            iterCt = iterCt + 1;
        end

        y(:,i+1) = yNew;
    end

end