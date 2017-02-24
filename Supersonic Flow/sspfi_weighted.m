function PHI = sspfi_weighted (XX, YY, yB, M0, w)
%% Supersonic Potential Flow - Implicit Solver
% Solves for flow potential in supersonic conditions using an implicit
% BVP solver. Includes a weighted correction to remove numerical noise.

% Get Beta and other important constants
beta = sqrt(abs(1 - M0^2));
dy = (YY(2,1)-YY(1,1));
dx = XX(1,2)-XX(1,1);

% Define Potential Field
PHI = XX;
[n_y, m_x] = size(XX);

for ii = 3:m_x
    a_vals = [0, (w/dy^2)*ones(1, n_y-2), (2*w/dy^2)];
    b_vals = (-(beta/dx)^2 - 2*w/dy^2)*ones(n_y);
    c_vals = [(2*w/dy^2), (w/dy^2)*ones(1, n_y-2), 0];
    d_vals = ((beta/dx)^2)*(-2*PHI(:,ii-1) + PHI(:,ii-2)) - ((1-w)/dy^2)*([0;PHI(1:end-2,ii-2);2*PHI(end-1,ii-2)] - 2.*PHI(:,ii-2) + [2*PHI(2,ii-2);PHI(3:end,ii-2);0]);
    
    % Check for boundary conditions
    if (XX(1,ii) > 0) && (XX(1,ii) <1)
        d_vals(1) = d_vals(1) + (1/(dx*dy))*(yB(ii) - yB(ii-2));
    end
    
    if (XX(1,ii-2) > 0) && (XX(1,ii-2) <1)
        d_vals(1) = d_vals(1) + (1/(dx*dy))*(yB(ii-2) - yB(ii-2-2));
    end
    
    PHI(:,ii) = thomas3(a_vals, b_vals, c_vals, d_vals);
end

end