function PHI = sspf_implicit (XX, YY, yB, M0)
%% Supersonic Potential Flow - Implicit Solver
% Solves for flow potential in supersonic conditions using an implicit
% BVP solver.

% Get Beta and other important constants
beta = sqrt(abs(1 - M0^2));
dy = (YY(2,1)-YY(1,1));
dx = XX(1,2)-XX(1,1);

% Define Potential Field
PHI = XX;
[n_y, m_x] = size(XX);

for ii = 3:m_x
    a_vals = [0, (1/dy^2)*ones(1, n_y-2), (2/dy^2)];
    b_vals = (-(beta/dx)^2 - 2/dy^2)*ones(n_y);
    c_vals = [(2/dy^2), (1/dy^2)*ones(1, n_y-2), 0];
    d_vals = ((beta/dx)^2)*(-2*PHI(:,ii-1) + PHI(:,ii-2));
    
    % Check for boundary conditions
    if (XX(1,ii) > 0) && (XX(1,ii) <1)
        d_vals(1) = d_vals(1) + (1/(dx*dy))*(yB(ii) - yB(ii-2));
    end
    
    PHI(:,ii) = thomas3(a_vals, b_vals, c_vals, d_vals);
end

end