function [PHI, PHI_R, PHI_T, RHO, res] = flowSolve_comp(grid, flow, times)

% Initilize Grid
dr = grid.dr;
dT = grid.dT;
dt = times.dt;
RR = grid.RR;
TT = grid.TT;
alpha = times.alpha;

% initalize grid count
n_t = (times.stop - times.start)/dt;
[n_r, n_T] = size(grid.RR);

% Initialize Values for flow field
PHI = zeros(n_r, n_T, n_t);
PHI(:,:,1) = RR.*cos(TT);
PHI(:,:,2) = PHI(:,:,1);

% Initialize flow field derivatives
RHO = zeros(size(PHI));
PHI_R = zeros(size(PHI));
PHI_T = zeros(size(PHI));
[PHI_R(:,:,2), PHI_T(:,:,2), RHO(:,:,2)] = update_derivs(PHI, RR, TT, flow.M0, flow.gamma);


% Initialize error checks
res = zeros(1,n_t-2);

% Begin iteration
for n = 3:n_t % loop through time
    
    for i = 1:n_T % loop through theta        
        for j = 1:(n_r-1) % loop through radius
            % 5 different control indexes for i and j directions            
            i0 = i;
            j0 = j;
            % Boundary Conditions for Theta
            if i == 1
                i1 = i+1;
                i_1 = i1;
            elseif i == n_T
                i0 = i;
                i1 = i-1;
                i_1 = i-1;
            else
                i1 = i+1;
                i_1 = i-1;
            end

            % Boundary Conditions for Radius
            if j == 1
                j1 = j+1;
                j_1 = j1;
            elseif j == n_r
                j1 = j+1;
                j_1 = j-1;
            else
                j1 = j+1;
                j_1 = j-1;
            end
            
            % Density averages
            RHO_i1 = 0.5*(RHO(j0,i1) + RHO(j0,i0));
            RHO_i0 = 0.5*(RHO(j0,i_1) + RHO(j0,i0));
            RHO_j1 = 0.5*(RHO(j1, i0) + RHO(j0,i0));
            RHO_j0 = 0.5*(RHO(j_1, i0) + RHO(j0,i0));
            
            % Second Order Derivatives
            PHI_TT = (1/RR(j0,i0)^2)*(1/dT)*(RHO_i1*(PHI(j0,i1,n-1)-PHI(j0,i0,n-1))/dT - RHO_i0*(PHI(j0,i0,n-1)-PHI(j0,i_1,n-1))/dT);
            PHI_RR = (1/RR(j0,i0))*((0.5*(RR(j1,i0)+RR(j0,i0))*(PHI(j1,i0,n-1)-PHI(j0,i0,n-1))/dr)*RHO_j1 - RHO_j0*(0.5*(RR(j0,i0)+RR(j_1,i0))*(PHI(j0,i0,n-1)-PHI(j_1,i0,n-1))/dr))/dr;

            % Apply Conditions for Time
            PHI(j0,i0,n) = ((0.5*alpha/dt - 1/dt^2)*PHI(j0,i0,n-2) + 2*PHI(j0,i0,n-1)/dt^2 + PHI_RR + PHI_TT)/(1/dt^2 + 0.5*alpha/dt);
        end
        PHI(n_r,i,n) = PHI(n_r,i,n-1); % set the boundary condition
    end
    
    % Update derivatives
    [PHI_R(:,:,n), PHI_T(:,:,n), RHO(:,:,n)] = update_derivs(PHI, RR, TT, flow.M0, flow.gamma);
    
    % Run Residual Check
    difference = abs(PHI(:,:,n) - PHI(:,:,n-1));
    res(n-2) = max(difference(:));
end

function [phi_r, phi_th, rho] = update_derivs(phi, rr, tt, M0, gam)
    
if ~exist('n', 'var')
    n = 2;
end

% initialize temp vals for PHI_R, PHI_T, and rho
phi_r = zeros(size(phi(:,:,n)));
phi_th = zeros(size(phi(:,:,n)));
rho = zeros(size(phi(:,:,n)));
% phi_t = zeros(size(rho));

    % Radial derivative
    for jj = 2:(n_r-1)
       phi_r(jj,:) =  (phi(jj+1,:,n)-phi(jj-1,:,n))./(rr(jj+1,:) - rr(jj-1,:));
    end

    phi_r(1,:) = zeros(1,n_T);
    phi_r(end,:) = (phi(end,:,n)-phi(end-1,:,n))./(rr(end,:) - rr(end-1,:));

    % Angular derivative
    for ii = 1:n_T
        ii0 = ii;
        if ii == 1
            ii1 = ii+1;
            ii_1 = n_T;
        elseif ii == n_T
            ii1 = 1;
            ii_1 = ii-1;
        else
            ii1 = ii+1;
            ii_1 = ii-1;
        end 

        phi_th(:,ii0) =  (1./rr(:,ii0)).*(phi(:,ii1,n)-phi(:,ii_1,n))./(tt(:,ii1) - tt(:,ii_1));        
    end

    % Time Derivative
    phi_t = (phi(:,:,n) - phi(:,:,n-1))/dt;

    % Local speed of sound
    a2_ij = (2/(M0^2) - (gam - 1).*(phi_r.^2 + (phi_th./rr).^2 - 1));

    % Local Mach number & artificial viscosity
    M_ij = (((phi_th./rr).^2 + phi_r.^2)./a2_ij).^2;
    eps = max(zeros(size(M_ij)), 1 - 0.9./(M_ij.^2));

    % Local density calculations
    rho_bar = (1 - 0.5.*(gam-1).*M0^2.*(phi_r.^2 + (phi_th./rr).^2 + 2*phi_t - 1)).^(1/(gam-1));

    for ii = 1:n_T % loop through Thetas
        ii0 = ii;
        if ii == 1
            ii_1 = n_T;
        else
            ii_1 = ii-1;
        end
        for jj = 2:n_r % Loop through Radii
            rho(jj,ii0) = rho_bar(jj, ii0) - eps(jj, ii0)*(phi_th(jj, ii0)*(rho_bar(jj,ii0) - rho_bar(jj,ii_1)) + phi_r(jj,ii0)*(rho_bar(jj,ii0) - rho_bar(jj-1, ii0)));
        end

        rho(1,ii0) = rho_bar(1,ii0) - eps(1,ii0) * (phi_th(1,ii0) * (rho_bar(1,ii0) - rho_bar(1, ii_1))); % all PHI_R's at surface are zero
    end
end

end

