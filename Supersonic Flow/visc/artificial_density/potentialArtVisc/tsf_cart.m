function [OUT] = tsf_cart(GR, FL, BC, CT)
%% INPUTS

% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
    % Types of inputs needed
    % XX, YY
    % dx, dy
    
% FL - fluid parameters
    % M0 - freestream mach number
    % gam - gamma of the flow

% CT - simulation control, including tolerances, viscous factor gain, etc.
    % Types of inputs needed
    % tol - tolerances
    % v_coeff - gain for artificial viscosity
    % t_rho - add time correction to density
    % alpha - damping for three-level iteration scheme
    % CT.dt - spacing for three-level scheme
    
% CT.tol = 1e-5;
% CT.v_coeff = 1.0;
% CT.t_rho = 1.0;


% BC - boundary conditions, including at wall conditions and far-field
% conditions
% % Boundary Conditions
% BC.Vy_II = zeros(size(GR.YY(end,:)));
% BC.Vx_II = ones(size(GR.XX(end,:)));
% BC.Vx_I = ones(size(GR.XX(:,1)));
% BC.PHI_II = GR.XX(end,:);
% BC.PHI_I = GR.XX(:,1);


%% OUTPUTS

% OUT - contains the resulting flow potential field, as well as the
% velocity, residuals, and indices of residual locations
    % PHI
    % U_n - 
    % M_ij - local mach number
    % a_l - local speed of sound

%% INITIALIZE

PHI = zeros([size(GR.XX), 2]);
PHI(:,:,1) = GR.XX;
PHI(:,:,2) = PHI(:,:,1);
PHI_new = PHI(:,:,2);
res = 1;
visc_on = 0;
ind = [];

mid = [-1/(GR.dy^2), (-2/GR.dy^2).*ones(1, length(GR.y_vals)-2), 0];
top = (1/GR.dy^2).*ones(1, length(GR.y_vals)-1);
bot = [(1/GR.dy^2).*ones(1, length(GR.y_vals)-2), 0];

sys_yy = diag(mid) + diag(top, 1) + diag(bot, -1);

%% START SOLUTION

while (res(end) > CT.tol)|| (length(res) < 100) % iterate through time

    % Reorganize three-level scheme
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI_new;

    % Apply boundary conditions
    PHI_new(:,1) = BC.PHI_I; % inlet condition

    % Calculate local viscosity at the nodes
    U_n = ones(size(GR.XX)); % PHI_X... only need velocity in X-dir
    for i = 2:(size(GR.XX, 2)-1) % loop through x-vals
        U_n(:,i) = (PHI(:,i+1,2) - PHI(:,i-1,2))./(2.*GR.dx);
    end
    U_n(:,1) = BC.Vx_I;
    if FL.M0 > 1
       U_n(:,end) =  U_n(:,end-1);
    end

    a2_ij = (1./FL.M0.^2)-0.5*(FL.gam-1).*(U_n.^2 - 1);
    eps_ij = max(zeros(size(U_n)), 1 - 1.0.*(a2_ij./(U_n.^2)));

    if any(any(eps_ij ~= 0)) && (visc_on ==0) % viscosity not being used and then turned on
        visc_on = 1;
        fprintf('Viscosity model activated! Iteration: %i\n', length(res));        
    elseif (visc_on == 1) && all(all(eps_ij == 0))
        visc_on = 0;
        fprintf('Viscosity model turned off! Iteration: %i\n', length(res));
    end
    
    % Check CFL condition
    V_n = diff(PHI(:,:,2),1,1)./GR.dy;
    CFL_i = (max(abs(U_n(:)))./GR.dx + max(max(abs(V_n(:))), max(BC.dyBdx))./GR.dy)*CT.dt;
    
    if CFL_i >= 1.0
       fprintf('CFL condition not met! Decreasing time steps!\n');
       CT.dt = CT.dt*0.8 / CFL_i;
       fprintf('New time step:%0.5f\n', CT.dt);
    end

    % Initialize Density calculations using between grid in X-dir
    u_avg = diff(PHI(:,:,2),1,2)./GR.dx;
    if CT.t_rho == 0
        rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(u_avg.^2 - 1)).^(1./(FL.gam-1));
    
    elseif CT.t_rho == 1 % add time correction
        phi_t = (PHI(:,:,2) - PHI(:,:,1))./CT.dt;
        phi_t_avg = 0.5.*(phi_t(:,1:(end-1)) + phi_t(:,2:end)); % average phi_t for three level scheme
        rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(u_avg.^2 + 2.*phi_t_avg - 1)).^(1./(FL.gam-1)); % do I need the time correction?
    end

    if any(imag(rho_ij(:)) ~= 0)
%         dbstop;
        test = imag(rho_ij);
        [y_ind, x_ind] = find(test ~= 0);
        fprintf('Non-real result for density!\n'); 
        fprintf('Check Nodes...\n');
        fprintf('Row %i, Col %i \n', y_ind, x_ind);
        fprintf('Limitting Density to Zero! \n');
        rho_ij(y_ind, x_ind) = 0;
    end

    % Calculate viscosity corrections
    rho_ds = zeros(size(rho_ij));
    for i = 1:(size(rho_ij,2)) % loop through x-dir
        if i ==1 % Dirichlet BC
            rho_ds(:,i) = 0.5*(rho_ij(:,i+1)-ones(size(rho_ij(:,i)))) - 0.5.*sign(u_avg(:,i)).*(rho_ij(:,i+1)-2.*rho_ij(:,i) + ones(size(rho_ij(:,i))));
        elseif i == size(rho_ij,2)
            rho_ds(:,i) = 0.5*(-rho_ij(:,i-1)+ones(size(rho_ij(:,i)))) - 0.5.*sign(u_avg(:,i)).*(rho_ij(:,i-1)-2.*rho_ij(:,i) + ones(size(rho_ij(:,i))));
        else
            rho_ds(:,i) = 0.5*(rho_ij(:,i+1)-rho_ij(:,i-1)) - 0.5.*sign(u_avg(:,i)).*(rho_ij(:,i+1)-2.*rho_ij(:,i) + rho_ij(:,i-1));
        end
    end

    if (FL.M0 == 0) && any(any(rho_ij ~= 1))
        fprintf('Incompressible Case not satisfied!\n');
    end

    if FL.M0 < 1 % subsonic, use Elliptic Solver
        for i = 2:(size(PHI,2)-1) % march along x-dir to calculate one step ahead, do not need to calculate initial inlet
            % Get Density
            if i == size(PHI,2)
                RHO_i1 = ones(size(rho_ij(:,i)));
            else
                RHO_i1 = rho_ij(:,i) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i); % use local viscosity at current node
            end
            RHO_i_1 = rho_ij(:,i-1) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-1);

            PHI_XX = (RHO_i1.*(PHI(:,i+1,2) - PHI(:,i,2))./GR.dx - RHO_i_1.*(PHI(:,i,2) - PHI(:,i-1,2))./GR.dx)./GR.dx;

            PHI_YY = sys_yy*PHI(:,i,2) - (1/GR.dy).*[BC.dyBdx(i); zeros(length(GR.y_vals)-1, 1)];

            PHI_new(:,i) = (PHI_XX + PHI_YY + 2/(CT.dt^2)*PHI(:,i,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,i,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
        end

        PHI_new(end,:) = BC.PHI_II;

    elseif FL.M0 > 1 % supersonic, use Hyperbolic Solver?

        for i = 3:(size(PHI,2)) % march along x-dir to calculate one step ahead, do not need to calculate initial inlet
            if i == size(PHI,2)
                RHO_i1 = rho_ij(:,i-1) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-1); % use local viscosity at current node
                RHO_i_1 = rho_ij(:,i-2) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-2);

                PHI_XX = (RHO_i1.*(PHI(:,i,2) - PHI(:,i-1,2))./GR.dx - RHO_i_1.*(PHI(:,i-1,2) - PHI(:,i-2,2))./GR.dx)./GR.dx;
                
%                 
%                 PHI_YY = sys_yy*PHI(:,i,2) - (1/GR.dy).*[BC.dyBdx(i); zeros(length(GR.y_vals)-1, 1)];
                PHI_YY = [((PHI(2:end,i,2) - PHI(1:end-1,i,2))./GR.dy - [BC.dyBdx(i); (PHI(2:end-1,i,2) - PHI(1:end-2,i,2))./GR.dy])./GR.dy; 0];

                PHI_new(:,i) = (PHI_XX + PHI_YY + 2/(CT.dt^2)*PHI(:,i,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,i,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
%                 PHI_new(:,i) = PHI_new(:,i-1);
            else % start from line three

                RHO_i1 = rho_ij(:,i) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i); % use local viscosity at current node
                RHO_i_1 = rho_ij(:,i-1) - CT.v_coeff.*eps_ij(:,i).*rho_ds(:,i-1);

                PHI_XX = (RHO_i1.*(PHI(:,i+1,2) - PHI(:,i,2))./GR.dx - RHO_i_1.*(PHI(:,i,2) - PHI(:,i-1,2))./GR.dx)./GR.dx;

%                 PHI_YY = sys_yy*PHI(:,i,2) - (1/GR.dy).*[BC.dyBdx(i); zeros(length(GR.y_vals)-1, 1)];
                PHI_YY = [((PHI(2:end,i,2) - PHI(1:end-1,i,2))./GR.dy - [BC.dyBdx(i); (PHI(2:end-1,i,2) - PHI(1:end-2,i,2))./GR.dy])./GR.dy; 0];

                PHI_new(:,i) = (PHI_XX + PHI_YY + 2/(CT.dt^2)*PHI(:,i,2) - (1/CT.dt^2 - 0.5*CT.alpha/CT.dt)*PHI(:,i,1))/(1/CT.dt^2 + 0.5*CT.alpha/CT.dt);
            end
        end
        PHI_new(end,:) = BC.PHI_II;
    end

    % Run Residual Check
    difference = abs(PHI_new(:,:) - PHI(:,:,2));
    [res(end+1), ind(end+1)] = max(difference(:));

    % Initialize next time step
    if (length(res) > 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
        figure(1);semilogy(1:length(res), res);
    end

end

%% POST PROCESS

U_n = ones(size(GR.XX)); % PHI_X... only need velocity in X-dir
PHI(:,:,2) = PHI_new;

for i = 2:(size(GR.XX, 2)-1) % loop through x-vals
    U_n(:,i) = (PHI(:,i+1,end) - PHI(:,i-1,end))./(2.*GR.dx);
    % determine velocity at the far-field end?
end
U_n(:,1) = BC.Vx_I;

if FL.M0 > 1
   U_n(:,end) =  U_n(:,end-1);
end

% Calculate uncorrected density 
%     rho_ij = (1 - 0.5*(FL.gam-1).*M  0.^2.*(u_avg.^2 - 1)).^(1./(FL.gam-1));
rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(U_n.^2 - 1)).^(1./(FL.gam-1));
a2_avg = (1./FL.M0.^2)-0.5*(FL.gam-1).*(U_n.^2 - 1);
M_ij = U_n./sqrt(a2_avg);

% Assign to output
OUT.PHI = PHI;
OUT.U_n = U_n;
OUT.a_l = sqrt(a2_avg);
OUT.rho_ij = rho_ij;
OUT.M_ij = M_ij;
OUT.res = res;
end