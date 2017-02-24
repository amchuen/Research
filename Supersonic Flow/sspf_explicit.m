function PHI = sspf_explicit (XX, YY, yB, M0)
%% Supersonic Potential Flow - Implicit Solver
% Solves for flow potential in supersonic conditions using an implicit
% BVP solver.

% Get Beta and other important constants
beta = sqrt(abs(1 - M0^2));
dy = (YY(2,1)-YY(1,1))*beta;
dx = XX(1,2)-XX(1,1);

% Define Potential Field
PHI = XX;
[n_y, m_x] = size(XX);

for ii = 2:m_x-1
    for jj = 1:n_y-1
        if jj == 1
            if (XX(jj,ii) > 0) && (XX(jj,ii) < 1) % body condition
                PHI(jj,ii+1) = (dx^2 / dy^2) * (2*PHI(jj+1, ii) - (dy/dx)*(yB(ii+1) - yB(ii-1)) - 2 * PHI(jj,ii)) + 2*PHI(jj,ii) - PHI(jj, ii-1);
            else
                PHI(jj,ii+1) = (dx^2 / dy^2) * (2*PHI(jj+1,ii) - 2 * PHI(jj,ii)) + 2*PHI(jj,ii) - PHI(jj, ii-1);
            end
        else
            PHI(jj, ii+1) = (dx^2 / dy^2) * (PHI(jj+1, ii) - 2 * PHI(jj,ii) + PHI(jj-1,ii)) + 2*PHI(jj,ii) - PHI(jj, ii-1);
        end
    end
    
end


end