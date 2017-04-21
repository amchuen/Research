function [PHI, PHI_Y, PHI_X, RHO, res, ind, x_res, y_res] = flowSolve_cartesian(grid, flow, times)

% Initilize Grid
dy = grid.dy;
dx = grid.dx;
dt = times.dt;
YY = grid.YY;
XX = grid.XX;
alpha = times.alpha;

% initalize grid count
n_t = (times.stop - times.start)/dt;
[n_y, n_x] = size(grid.XX);

% Initialize Values for flow field
PHI = zeros(n_y, n_x, n_t);
PHI(:,:,1) = flow.PHI_INIT;
PHI(:,:,2) = PHI(:,:,1);

% Boundary Conditions for field at time n
PHI_BC = grid.XX; %YY(n_y,:).*cos(XX(n_y,:)) + cos(XX(n_y,:))./YY(n_y,:);
% PHI_cyl = RR(1,:).*cos(TT(1,:)) + cos(TT(1,:))./RR(1,:);

% Initialize flow field derivatives
RHO = zeros(size(PHI));
PHI_Y = zeros(size(PHI));
PHI_X = zeros(size(PHI));
[PHI_Y(:,:,2), PHI_X(:,:,2), RHO(:,:,2)] = update_rho(PHI, YY, flow.M0, flow.gamma);
e_visc = -0.001;


% Initialize error checks
res = zeros(1,n_t-2);
x_res = res;
y_res = res;
ind = res;

% Begin iteration
for n = 3:n_t % loop through time
%     if (mean2(RHO(:,:,n-1)) ~= 1)
%        fprintf('Warning: density is not unity!\n'); 
%     end
        
    for i = 1:n_x % loop through theta        
        for j = 1:(n_y-1) % loop through radius
            % 5 different control indexes for i and j directions            
            i0 = i;
            j0 = j;
            % Boundary Conditions for Theta
            if i == 1
                i1 = i+1;
                i_1 = i1;
            elseif i == n_x
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
            else
                j1 = j+1;
                j_1 = j-1;
            end
            
            % Density averages
            RHO_i1 = 0.5*(RHO(j0,i1, n-1) + RHO(j0,i0, n-1));
            RHO_i0 = 0.5*(RHO(j0,i_1,n-1) + RHO(j0,i0,n-1));
            RHO_j1 = 0.5*(RHO(j1, i0,n-1) + RHO(j0,i0,n-1));
            RHO_j0 = 0.5*(RHO(j_1, i0,n-1) + RHO(j0,i0,n-1));
            
            if (flow.M0 == 0) && ((RHO_i1 ~= 1) || (RHO_i0 ~= 1) || (RHO_j1 ~= 1) || (RHO_j0 ~= 1) )
               fprintf('WARNING: DENISTY IS NOT UNITY AT ITER %i, %i, %i\n', n, i, j);  
            end
            
            % Divergence of Field
            PHI_XX = (1/dx)*(RHO_i1*(PHI(j0,i1,n-1)-PHI(j0,i0,n-1))/dx - RHO_i0*(PHI(j0,i0,n-1)-PHI(j0,i_1,n-1))/dx);
            PHI_YY = (1/dy)*(RHO_j1*(PHI(j1,i0,n-1) - PHI(j0,i0,n-1))/dy - RHO_j0*(PHI(j0,i0,n-1) - PHI(j_1,i0,n-1))/dy);
            
            % Divergence of Density
            RHO_XX = (1/dx^2)*(RHO(j0,i1,n-1) - 2*RHO(j0,i0,n-1)+RHO(j0,i_1,n-1));
            RHO_YY = ((RHO(j1,i0,n-1) - RHO(j0,i0,n-1))/dy - (RHO(j0,i0,n-1) - RHO(j_1,i0,n-1))/dy)/dy;
            
            if (j == 1) && (PHI(j1,i0) ~= PHI(j_1,i0))
               fprintf('Neumann Condition not applied for radius!\n'); 
            end
            
            if ((i == 1)||(i == n_x)) && (PHI(j0,i1) ~= PHI(j0,i_1))
               fprintf('Neumann Condition not applied for THETA!\n'); 
            end

            % Apply Conditions for Time ...add viscosity term from divergence of density
            PHI(j0,i0,n) = ((0.5*alpha/dt - 1/(dt^2))*PHI(j0,i0,n-2) + 2/(dt^2)*PHI(j0,i0,n-1) + PHI_YY + PHI_XX + e_visc.*(RHO_YY + RHO_XX))/(1/dt^2 + 0.5*alpha/dt);
        end
    end
    PHI(n_y,:,n) = PHI_BC; %PHI(n_r,i,n-1); % set the farfield boundary condition
    
    % Update derivatives
    [PHI_Y(:,:,n), PHI_X(:,:,n), RHO(:,:,n)] = update_rho(PHI, YY, flow.M0, flow.gamma);
    
    % Run Residual Check
    difference = abs(PHI(:,:,n) - PHI(:,:,n-1));
    [res(n-2), ind(n-2)] = max(difference(:));
    [yy, xx] = ind2sub(size(difference), ind(n-2));
    y_res(n-2) = yy;
    x_res(n-2) = xx;
end

function [phi_r, phi_th, rho] = update_rho(phi, rr, M0, gam)
    
    if ~exist('n', 'var')
        n = 2;
    end

    % initialize temp vals for PHI_R, PHI_T, and rho
    phi_r = zeros(size(phi(:,:,n)));
    phi_th = zeros(size(phi(:,:,n)));
    rho = ones(size(phi(:,:,n)));
    % phi_t = zeros(size(rho));

    % Radial derivative
    for jj = 2:(n_y-1)
       phi_r(jj,:) =  (phi(jj+1,:,n)-phi(jj-1,:,n))./(2*dy); % central finite difference
    end

    phi_r(1,:) = zeros(size(phi_r(1,:))); % surface boundary condition
    phi_r(end,:) = (phi(end,:,n)-phi(end-1,:,n))./(dy); % backwards

    % Angular derivative
    for ii = 2:(n_x-1)
        ii0 = ii;
        ii1 = ii+1;
        ii_1 = ii-1;
        % Central differencing scheme
        phi_th(:,ii0) =  (1./rr(:,ii0)).*(phi(:,ii1,n)-phi(:,ii_1,n))./dx;        
    end

    phi_th(:,1) = zeros(size(phi_th(:,1))); % symmetric boundary condition
    phi_th(:,n_x) = zeros(size(phi_th(:,n_x))); % symmetric boundary condition

    % Time Derivative
    phi_t = (phi(:,:,n) - phi(:,:,n-1))/dt;

    % Local speed of sound
    a2_ij = (2/(M0^2) - (gam - 1).*(phi_r.^2 + (phi_th./rr).^2 - 1));

    % Local Mach number & artificial viscosity
    M_ij = (((phi_th./rr).^2 + phi_r.^2)./a2_ij).^2;
    eps = max(zeros(size(M_ij)), 1 - 0.9./(M_ij.^2));

    % Local density calculations
    rho_bar = (1 - 0.5.*(gam-1).*M0^2.*(phi_r.^2 + (phi_th./rr).^2 + 2*phi_t - 1)).^(1/(gam-1));

    for ii = 1:n_x % loop through Thetas
        ii0 = ii;
        if ii == 1
            ii_1 = n_x;
        else
            ii_1 = ii-1;
        end
        for jj = 2:(n_y-1) % Loop through Radii
            rho(jj,ii0) = rho_bar(jj, ii0) - eps(jj, ii0)*(phi_th(jj, ii0)*(rho_bar(jj,ii0) - rho_bar(jj,ii_1)) + phi_r(jj,ii0)*(rho_bar(jj,ii0) - rho_bar(jj-1, ii0)));
            if rho(jj,ii0) < 0
                rho(jj,ii0) = 0;
            end
        end

        rho(1,ii0) = rho_bar(1,ii0) - eps(1,ii0) * (phi_th(1,ii0) * (rho_bar(1,ii0) - rho_bar(1, ii_1))); % all PHI_R's at surface are zero
        if rho(1,ii0) < 0
            rho(1,ii0) = 0;
        end
        
        rho(n_y, ii0) = 1;
    end
end

end

