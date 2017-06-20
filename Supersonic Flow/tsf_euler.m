function [OUT] = tsf_euler(GR, FL, BC, CT)
%% INPUTS

% GR - grid information, such as the meshfield, grid spacing (dx, dy, etc.)
    % Types of inputs needed
    % RR, TT
    % dr, dT
    
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
    % eps - viscosity coefficient?
    
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

PHI = zeros([size(GR.RR), 3]);
PHI(:,:,1) = GR.RR;
PHI(:,:,2) = PHI(:,:,1);
PHI(:,:,3) = PHI(:,:,2);
RHO = ones([size(GR.RR), 3]);
res = 1;
visc_on = 0;
ind = [];

% mid = [-1/(GR.dy^2), (-2/GR.dy^2).*ones(1, length(GR.y_vals)-2), 0];
% top = (1/GR.dy^2).*ones(1, length(GR.y_vals)-1);
% bot = [(1/GR.dy^2).*ones(1, length(GR.y_vals)-2), 0];
% 
% sys_yy = diag(mid) + diag(top, 1) + diag(bot, -1);

%% START SOLUTION

while (res(end) > CT.tol)|| (length(res) < 100) % iterate through time
        
    % Reorganize three-level scheme
    PHI(:,:,1) = PHI(:,:,2);
    PHI(:,:,2) = PHI(:,:,3);
    
    RHO(:,:,1) = RHO(:,:,2);
    RHO(:,:,2) = RHO(:,:,3);

    % Apply boundary conditions
%     PHI(:,1,3) = BC.PHI_I; % inlet condition
    
%     Calculate Velocities
    phiT = [zeros(size(GR.RR(:,1))),... % Neumann condition at LE/TE of cyl b/c of symmetry condition
            (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*GR.dT),... % Central Difference
            zeros(size(GR.RR(:,end)))]; % PHI_X
        
    phiR  = [zeros(size(PHI(1,:, 2)));... % Neumann condition at cylinder surface
			(PHI(3:end,:,2) - PHI(1:(end-2),:,2))./(2*GR.dr);...
			BC.Vr_II]; % PHI_Y, central difference
        
    q2_ij = (phiT ./ GR.RR).^2 + phiR.^2;
%     a2_ij = (1./FL.M0.^2)-0.5*(FL.gam-1).*(q2_ij - 1);
% 	eps_ij = max(zeros(size(a2_ij)), 1 - 0.99.*(a2_ij./q2_ij));

%     if any(eps_ij(:) ~= 0) && (visc_on ==0) % viscosity not being used and then turned on
%         visc_on = 1;
%         fprintf('Viscosity model activated! Iteration: %i\n', length(res)+1);        
%     elseif (visc_on == 1) && (all(eps_ij(:) == 0))
%         visc_on = 0;
%         fprintf('Viscosity model turned off! Iteration: %i\n', length(res)+1);
%     end
    
    for i = 1:(size(PHI,2)) % march along theta direction
        % Density Calculation
        % Calculate all values except north boundary, which will be set as
        % 1 for non-dimensionalized far-field density
        if (i == 1)
            i1 = i+1;
            i_1 = i1;
        elseif(i == size(PHI,2))
            i_1 = i-1;
            i1 = i_1;
        else
            i_1 = i-1;
            i1 = i+1;
        end
        
        rho_laplace = (1./GR.RR(1:(end-1),i)).^2 .* ((RHO(1:(end-1),i_1,2) - 2.*RHO(1:(end-1),i,2) + RHO(1:(end-1),i1,2))./(GR.dT^2)) + ...
                       (1./GR.RR(1:(end-1),i)).*((0.5.*(GR.RR(2:(end),i)+GR.RR(1:(end-1),i)).*diff(RHO(:,i,2))... ([RHO(2:end,i,2) - RHO(1:end-1,i,2)])
                       - 0.5.*([2*GR.RR(1,i); GR.RR(1:(end-2),i)+GR.RR(2:(end-1),i)]).*([0; diff(RHO(1:end-1,i,2))...RHO(2:end-1,i,2) - RHO(1:end-2,i,2)
                       ]))./(GR.dr^2));
        phi_laplace = (1./GR.RR(1:(end-1),i)).^2 .* ((PHI(1:(end-1),i_1,2) - 2.*PHI(1:(end-1),i,2) + PHI(1:(end-1),i1,2))./(GR.dT^2)) + ...
                       (1./GR.RR(1:(end-1),i)).*((0.5.*(GR.RR(2:(end),i)+GR.RR(1:(end-1),i)).*diff(PHI(:,i,2))... ([RHO(2:end,i,2) - RHO(1:end-1,i,2)])
                       - 0.5.*([2*GR.RR(1,i); GR.RR(1:(end-2),i)+GR.RR(2:(end-1),i)]).*([0; diff(PHI(1:end-1,i,2))...RHO(2:end-1,i,2) - RHO(1:end-2,i,2)
                       ]))./(GR.dr^2));
        
        RHO(:,i,3) = [RHO(1:(end-1),i,2) + 2.*CT.dt.*(CT.v_coeff.*rho_laplace - ...
                        (1./(GR.RR(1:(end-1),i).^2)).*((0.5.*(RHO(1:(end-1),i1,2) + RHO(1:(end-1),i,2)).*(PHI(1:(end-1),i1,2) - PHI(1:(end-1),i,2))...
                                                    - 0.5.*(RHO(1:(end-1),i,2)+RHO(1:(end-1),i_1,2)).*(PHI(1:(end-1),i,2)-PHI(1:(end-1),i_1,2)))./(GR.dT^2)) - ...
                        (1./(GR.RR(1:(end-1),i))).*(0.5.*(RHO(2:(end),i,2)+RHO(1:(end-1),i,2)).*0.5.*(GR.RR(2:(end),i)+GR.RR(1:(end-1),i)).*(diff(PHI(:,i,2))) - ...
                                                    0.5.*([2*RHO(1,i,2); RHO(1:(end-2),i,2)+RHO(2:(end-1),i,2)]).*0.5.*([2*GR.RR(1,i); GR.RR(1:(end-2),i)+GR.RR(2:(end-1),i)]).*([0; diff(PHI(1:end-1,i,2))]))./(GR.dr^2));...
                       1];
                   
        PHI(:,i,3) = [PHI(1:(end-1),i,2) + 2.*CT.dt.*(CT.v_coeff.*phi_laplace...
                        - (RHO(1:(end-1),i,2).^(FL.gam-1))./((FL.gam-1)*FL.M0)...
                        -0.5*q2_ij(1:(end-1),i)...- 0.5.*(((1./GR.RR(1:(end-1),i)).*(()./)).^2 + ().^2)...
                        + 0.5 + (1 / ((FL.gam - 1)*FL.M0)));...
                        GR.RR(end,i).*cos(GR.TT(end,i))];
    
    end

    % Run Residual Check
    diff_phi = abs(PHI(:,:,3) - PHI(:,:,2));
    diff_rho = abs(RHO(:,:,3) - RHO(:,:,2));
    res(end+1)= max(max(diff_phi(:)), max(diff_rho(:)));

    % READOUT
    if (length(res) > 500) && (mod(length(res), 500) == 0)
        fprintf('Iteration Ct: %i\n', length(res));
        fprintf('Current Residual: %0.5e\n', res(end));
        figure(1);semilogy(1:length(res), res);
    end

end

%% POST PROCESS

% U_n = ones(size(GR.XX)); % PHI_X... only need velocity in X-dir
% PHI(:,:,2) = PHI_new;
% 
% for i = 2:(size(GR.XX, 2)-1) % loop through x-vals
%     U_n(:,i) = (PHI(:,i+1,end) - PHI(:,i-1,end))./(2.*GR.dx);
%     % determine velocity at the far-field end?
% end
% U_n(:,1) = BC.Vx_I;
% 
% if FL.M0 > 1
%    U_n(:,end) =  U_n(:,end-1);
% end

U_n = [zeros(size(GR.RR(:,1))),... % Neumann condition at LE/TE of cyl b/c of symmetry condition
            (PHI(:,3:end,end) - PHI(:,1:(end-2),end))./(2*GR.dT),... % Central Difference
            zeros(size(GR.RR(:,end)))]; % PHI_X
        
V_n = [zeros(size(PHI(1,:, 2)));... % Neumann condition at cylinder surface
			(PHI(3:end,:,2) - PHI(1:(end-2),:,2))./(2*GR.dr);...
			BC.Vr_II]; % PHI_Y, central difference

% Calculate uncorrected density 
%     rho_ij = (1 - 0.5*(FL.gam-1).*M  0.^2.*(u_avg.^2 - 1)).^(1./(FL.gam-1));
rho_ij = (1 - 0.5*(FL.gam-1).*FL.M0.^2.*(U_n.^2 - 1)).^(1./(FL.gam-1));
q2_ij = (phiT ./ GR.RR).^2 + phiR.^2;
a2_ij = (1./FL.M0.^2)-0.5*(FL.gam-1).*(q2_ij - 1);
M_ij = sqrt(q2_ij./a2_ij);

% Assign to output
OUT.PHI = PHI;
OUT.U_n = U_n;
OUT.V_n = V_n;
OUT.a_l = sqrt(a2_ij);
OUT.rho_ij = RHO(:,:,3);
OUT.M_ij = M_ij;
OUT.res = res;
end