clc;
clear;
close all;

%% SIM CONTROL PARAMS - GRID INITIALIZATION

% Step Sizes
dr = 0.1;
dT = 0.01*pi;
dt = 0.1*dr;

alpha = 0.05;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 10;

% Field Axis Values
r_vals = (r_cyl+0.5*dr):dr:(r_max);
T_vals = 0:dT:(pi);

% Time Inputs
% start = 0.0;
% stop = 100 * run.times.dt;
% dt = [0.001];

% Grid Setup
[TT, RR] = meshgrid(T_vals, r_vals);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% Fluid Params
gam = 1.4; % heat 
M0 = 0.01;
visc = 0.0;

%% FIELD VARIABLE INITIALIZATION

PHI = zeros([size(RR), 2]);
res = 1;
ind = 0;
iter = 1.0;
tol = 1e-3;

BC.dirichlet.Vr_II = cos(TT(end,:)).*(1 - (r_cyl^2)./(RR(end,:).^2));
BC.dirichlet.PHI_II = (RR(end,:) + (r_cyl^2)./(RR(end,:))).*cos(TT(end,:));

%% START SOLUTION

while res(end) > tol % iterate through time
    
    % Initialize next time step
    iter = iter + 1;
    PHI(:,:,iter+1) = zeros(size(RR));

    % Calculate Velocity as function of PHI(n)
    U_n = zeros(size(RR)); % PHI_T
    V_n = zeros(size(RR)); % PHI_R
    for i = 2:(length(T_vals)-1) % loop through theta
        U_n(:,i) = (PHI(:,i+1,iter) - PHI(:,i-1,iter))./(2.*dT .* RR(:, i));
    end

    for ii = 1:(length(r_vals)-1) % loop through radius
        V_n(ii,:) = (PHI(ii+1,:,iter) - PHI(ii,:,iter))./(dr);
    end
    V_n(end,:) = BC.dirichlet.Vr_II; % apply far-field potential flow condition
    
    % Get Local Viscosity
    q2_ij = U_n.^2 + V_n.^2;
    a2_ij = (1./M0.^2)-0.5*(gam-1).*(q2_ij - 1);
    M2_ij = q2_ij./a2_ij; % local Mach number
    eps_ij = max(zeros(size(M2_ij)), 1 - 1.0./(M2_ij));
    
    % Get Density
    phi_t = (PHI(:,:,iter) - PHI(:,:,iter-1))./dt;
    rho_ij = (1 - 0.5*(gam-1).*M0.^2.*(q2_ij + 2*phi_t - 1)).^(1./(gam-1));
    del_rho_dT = zeros(size(rho_ij));
    for i = 1:(length(T_vals)) % loop through theta
        if i ==1 % Neumann BC
            i1 = i+1;
            i_1 = i1;
        elseif i == length(T_vals) % Neumann BC
            i1 = i-1;
            i_1 = i1;
        else
            i1 = i+1;
            i_1 = i-1;
        end
        del_rho_dT(:,i) = 0.5.*(rho_ij(:,i1)-rho_ij(:,i_1))-0.5.*(sign(U_n(:,i))).*(rho_ij(:,i1)-2.*rho_ij(:,i)+rho_ij(:,i_1));
    end
    
    del_rho_dr = zeros(size(rho_ij));
    for ii = 1:(length(r_vals)) % loop through radius
        if ii ==1 % Neumann BC
            del_rho_dr(ii,:) = 0.5*(rho_ij(ii+1,:) - rho_ij(ii,:))-0.5.*(sign(V_n(ii,:))).*(rho_ij(ii+1,:)-rho_ij(ii,:));
        elseif ii == length(r_vals) % Dirichlet BC
            del_rho_dr(ii,:) = zeros(size(rho_ij(end,:)));
        else
            del_rho_dr(ii,:) = 0.5*(rho_ij(ii+1,:) - rho_ij(ii,:))-0.5.*(sign(V_n(ii,:))).*(rho_ij(ii+1,:)-2.*rho_ij(ii,:)-rho_ij(ii-1,:));
        end
    end
    
    if all(max(q2_ij(1,:)) == 0)
        rho_ds = zeros(size(U_n));        
    else
        rho_ds = (U_n./sqrt(q2_ij)).*del_rho_dT + (V_n./sqrt(q2_ij)).*del_rho_dr;
    end
    RHO = rho_ij - eps_ij.*rho_ds;
    
    % Initialize Solution
    PHI_T = zeros(size(RR));
    PHI_R = zeros(size(RR));
    
    for i = 1:length(T_vals) % loop through theta        
        for j = 1:(length(r_vals)-1) % loop through radius
            % 5 different control indexes for i and j directions            
            i0 = i;
            j0 = j;
            % Boundary Conditions for Theta
            if i == 1
                i1 = i+1;
                i_1 = i1;
            elseif i == length(T_vals)
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
            RHO_i1 = 0.5*(RHO(j0,i1) + RHO(j0,i0));
            RHO_i0 = 0.5*(RHO(j0,i_1) + RHO(j0,i0));
            RHO_j1 = 0.5*(RHO(j1, i0) + RHO(j0,i0));
            RHO_j0 = 0.5*(RHO(j_1, i0) + RHO(j0,i0)); 
            
            % Radius Averages
            RR_i1 = 0.5*(RR(j0,i1) + RR(j0,i0));
            RR_i0 = 0.5*(RR(j0,i_1) + RR(j0,i0));
            RR_j1 = 0.5*(RR(j1, i0) + RR(j0,i0));
            RR_j0 = 0.5*(RR(j_1, i0) + RR(j0,i0));
            
            if j == 1
                PHI_RR = (RR_j1*RHO_j1*V_n(j0,i0))/(RR(j0,i0)*dr);
            else
                PHI_RR = (RR_j1*RHO_j1*V_n(j0,i0) - RR_j0*RHO_j0*V_n(j_1,i0))/(RR(j0,i0)*dr);
            end
            
            PHI_TT = (1/RR(j0,i0)^2)*(1/dT)*(RHO_i1*U_n(i0) - RHO_i0*U_n(i_1));
            
            PHI(j0,i0,iter+1) = ((0.5*alpha/dt - 1/(dt^2))*PHI(j0,i0,iter-1) + 2/(dt^2)*PHI(j0,i0,iter) + PHI_RR + PHI_TT)/(1/dt^2 + 0.5*alpha/dt);
            
        end
    end
    
    PHI(end,:,iter+1) = BC.dirichlet.PHI_II;
    
    % Run Residual Check
    difference = abs(PHI(:,:,iter+1) - PHI(:,:,iter));
    [res(iter), ind(iter)] = max(difference(:));
%     [yy, xx] = ind2sub(size(difference), ind(iter));
%     y_res(n-2) = yy;
%     x_res(n-2) = xx;
end

figure();
contourf(XX, YY, PHI(:,:,end));