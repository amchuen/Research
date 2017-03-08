clc;
clear;
close all;

%% SIM CONTROL PARAMS

% Step Sizes
dr = 0.1;
dT = 0.01*pi;
dt = 0.1*dr;
alpha = 5;

CFL = dt/dr + dt/dT;
disp(CFL);
if CFL >=1
    fprintf('CFL too high!\n');
end

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 5;

% Field Axis Values
r_vals = r_cyl:dr:(r_max+dr);
T_vals = 0:dT:(pi);


[TT, RR] = meshgrid(T_vals, r_vals);

XX = RR .* cos(TT);
YY = RR .* sin(TT);


[n_r, n_T] = size(RR);
n_t = 200;

% Initialize Values for flow field
PHI = zeros(n_r, n_T, n_t);
PHI(:,:,1) = RR.*cos(TT);
PHI(:,:,2) = PHI(:,:,1);

% Boundary Conditions for field at time n
PHI_BC = RR(n_r,:).*cos(TT(n_r,:)) + cos(TT(n_r,:))./RR(n_r,:);
PHI_cyl = RR(1,:).*cos(TT(1,:)) + cos(TT(1,:))./RR(1,:);


%% Algorithm

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
            else
                j1 = j+1;
                j_1 = j-1;
            end
            
            % Second Order Derivatives
            PHI_TT = (1/RR(j0,i0)^2)*(1/dT)*((PHI(j0,i1,n-1)-PHI(j0,i0,n-1))/dT - (PHI(j0,i0,n-1)-PHI(j0,i_1,n-1))/dT);
            PHI_RR = (1/RR(j0,i0))*((dr/2 + RR(j0,i0))*(PHI(j1,i0,n-1) - PHI(j0,i0,n-1))/dr - ((RR(j0,i0) - dr/2))*(PHI(j0,i0,n-1) - PHI(j_1,i0,n-1))/dr)/dr;
            
            if (j == 1) && (PHI(j1,i0) ~= PHI(j_1,i0))
               fprintf('Neumann Condition not applied for radius!\n'); 
            end
            
            if ((i == 1)||(i == n_T)) && (PHI(j0,i1) ~= PHI(j0,i_1))
               fprintf('Neumann Condition not applied for THETA!\n'); 
            end

            % Apply Conditions for Time
            PHI(j0,i0,n) = ((0.5*alpha/dt - 1/(dt^2))*PHI(j0,i0,n-2) + 2/(dt^2)*PHI(j0,i0,n-1) + PHI_RR + PHI_TT)/(1/dt^2 + 0.5*alpha/dt);
        end
    end
    PHI(n_r,:,n) = PHI_BC; %PHI(n_r,i,n-1); % set the farfield boundary condition
    
    % Check Derivatives
    difference = abs(PHI_cyl - PHI(1,:,n));
    
    % Update derivatives
%     [PHI_R(:,:,n), PHI_T(:,:,n), RHO(:,:,n)] = update_rho(PHI, RR, flow.M0, flow.gamma);
    
    % Run Residual Check
%     difference = abs(PHI(:,:,n) - PHI(:,:,n-1));
    res(n-2) = max(difference(:));
end

figure();
contourf(XX,YY,PHI(:,:,end), 50);

figure();
plot(1:(n_t - 2), res);