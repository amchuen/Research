clc;
clear;
close all;

%% FIELD SETUP

% Step Sizes
dr = 0.1;
dT = 0.01;
% dt = 0.1*dr*dT;

% Cylinder Dimensions
r_cyl = 1.0;
r_max = 3;

% Field Axis Values
r_vals = r_cyl:dr:(r_max+dr);
T_vals = 0:dT:(pi);


dts = [0.1, 0.05, 0.01]*dr*dT;
M_vals = [0];
gam = 1.4;

% Build Field
[TT, RR] = meshgrid(T_vals, r_vals);

XX = RR .* cos(TT);
YY = RR .* sin(TT);

% PHI = zeros(length(r_vals), length(T_vals), length(t_vals));

% Define the density equation


%% Boundary Conditions

% PHI(:,:,1) = RR.*cos(TT) + cos(TT)./RR;
% PHI(:,:,2) = PHI(:,:,1);

PHI_T1 = ones(length(M_vals), length(T_vals));
PHI_RHO1 = ones(size(PHI_T1));
res_time = {};

for tt = 1:length(dts)
    t_vals = 0:dts(tt):(100*dts(tt));
    dt = dts(tt);
    alpha = 10/(dt);
    res = ones(size(t_vals));
    PHI = zeros(length(r_vals), length(T_vals), length(t_vals));
    for mm = 1:length(M_vals)
        M0 = M_vals(1);
        PHI(:,:,1) = RR.*cos(TT) + cos(TT)./RR;
        PHI(:,:,2) = PHI(:,:,1);
        for n = 3:length(t_vals) % loop through time
            for i = 1:length(T_vals) % loop through theta        
                for j = 1:(length(r_vals)-1) % loop through radius
                    if j ~= (length(r_vals)-1) % only run this if not in the boundary condition
                        % 5 different control indexes for i and j directions            
                        i0 = i;
                        j0 = j;
                        % Boundary Conditions for Theta
                        if i == 1
                            i1 = i+1;
                            i2 = i+2;
                            i_1 = i1;
                            i_2 = i2;
                        elseif i == 2
                            i1 = i+1;
                            i2 = i+2;
                            i_1 = i-1;
                            i_2 = i;
                        elseif i == length(T_vals)
                            i0 = i;
                            i1 = i-1;
                            i2 = i-2;
                            i_1 = i-1;
                            i_2 = i-2;
                        elseif i == length(T_vals)-1
                            i1 = i+1;
                            i2 = i;
                            i_1 = i-1;
                            i_2 = i-2;
                        else
                            i1 = i+1;
                            i2 = i+2;
                            i_1 = i-1;
                            i_2 = i-2;
                        end

                        % Boundary Conditions for Radius
                        if j == 1
                            j1 = j+1;
                            j2 = j+2;
                            j_1 = j1;
                            j_2 = j2;
                            n1 = n-1; % use with j_1
                            n2 = n-1; % use with j_2
                        elseif j ==2
                            j1 = j+1;
                            j2 = j+2;
                            j_1 = j-1;
                            j_2 = j;
                            n1 = n-1; % use with j_1
                            n2 = n-1; % use with j_2
                        elseif j == length(r_vals)-3
                            j1 = j+1;
                            j2 = j+2;
                            j_1 = j-1;
                            j_2 = j-2;
                            n1 = n-1; % use with j_1
                            n2 = n-1; % use with j_2
                        elseif j == length(r_vals)-2
                            j1 = j+1;
                            j2 = j;
                            j_1 = j-1;
                            j_2 = j-2;
                            n1 = n-1; % use with j_1
                            n2 = n-1; % use with j_2
                        else
                            j1 = j+1;
                            j2 = j+2;
                            j_1 = j-1;
                            j_2 = j-2;
                            n1 = n-1; % use with j_1
                            n2 = n-1; % use with j_2
                        end

                        % RHO(j,i)
                        RHO_ij = (1 - 0.5*(gam-1)*M0^2*(((PHI(j1, i0, n-1) - PHI(j_1, i0, n1))/(2*dr))^2 + (1/RR(j0,i0))^2*((PHI(j0, i1, n-1)-PHI(j0,i_1,n-1))/(2*dT))^2 + 2*(PHI(j0,i0,n-1) - PHI(j0,i0,n-2))/dt - 1))^(1/(gam-1));

                        % RHO(j+1, i)
                        RHO_ij1 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j2, i0, n-1) - PHI(j0, i0, n-1))/(2*dr))^2 + (1/RR(j1,i0))^2*((PHI(j1, i1, n-1)-PHI(j1,i_1,n-1))/(2*dT))^2 + 2*(PHI(j1,i0,n-1) - PHI(j1,i0,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j-1, i)
                        RHO_ij0 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j0, i0, n-1) - PHI(j_2, i0, n2))/(2*dr))^2 + (1/RR(j_1,i0))^2*((PHI(j_1, i1, n1)-PHI(j_1,i_1,n1))/(2*dT))^2 + 2*(PHI(j_1,i0,n1) - PHI(j_1,i0,n1-1))/dt - 1))^(1/(gam-1));

                        % RHO(j,i+1)
                        RHO_i1j = (1 - 0.5*(gam-1)*M0^2*(((PHI(j1, i1, n-1) - PHI(j_1, i1, n1))/(2*dr))^2 + (1/RR(j0,i1))^2*((PHI(j0, i2, n-1)-PHI(j0,i0,n-1))/(2*dT))^2 + 2*(PHI(j0,i1,n-1) - PHI(j0,i1,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j,i-1)
                        RHO_i0j = (1 - 0.5*(gam-1)*M0^2*(((PHI(j1, i_1, n-1) - PHI(j_1, i_1, n1))/(2*dr))^2 + (1/RR(j0,i_1))^2*((PHI(j0, i0, n-1)-PHI(j0,i_2,n-1))/(2*dT))^2 + 2*(PHI(j0,i_1,n-1) - PHI(j0,i_1,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j+1,i+1)
                        RHO_i1j1 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j2, i1, n-1) - PHI(j0, i1, n-1))/(2*dr))^2 + (1/RR(j1,i1))^2*((PHI(j1, i2, n-1)-PHI(j1,i0,n-1))/(2*dT))^2 + 2*(PHI(j1,i1,n-1) - PHI(j1,i1,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j+1,i-1)
                        RHO_i0j1 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j2, i_1, n-1) - PHI(j0, i_1, n-1))/(2*dr))^2 + (1/RR(j1,i_1))^2*((PHI(j1, i0, n-1)-PHI(j1,i_2,n-1))/(2*dT))^2 + 2*(PHI(j1,i_1,n-1) - PHI(j1,i_1,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j-1,i+1)
                        RHO_i1j0 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j0, i1, n-1) - PHI(j_2, i1, n2))/(2*dr))^2 + (1/RR(j_1,i1))^2*((PHI(j_1, i2, n1)-PHI(j_1,i0,n1))/(2*dT))^2 + 2*(PHI(j_1,i1,n1) - PHI(j_1,i1,n1-1))/dt -1))^(1/(gam-1));

                        % RHO(j-1,i-1)
                        RHO_i0j0 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j0, i_1, n-1) - PHI(j_2, i_1, n2))/(2*dr))^2 + (1/RR(j_1,i_1))^2*((PHI(j_1, i0, n1)-PHI(j_1,i_2,n1))/(2*dT))^2 + 2*(PHI(j_1,i_1,n1) - PHI(j_1,i_1,n1-1))/dt -1))^(1/(gam-1));

                        % RHO(j, i-2) - (forward difference)
                        RHO_i_1j = (1 - 0.5*(gam-1)*M0^2*(((PHI(j1, i_2, n-1) - PHI(j_1, i_2, n1))/(2*dr))^2 + (1/RR(j0,i_2))^2*((PHI(j0, i_1, n-1)-PHI(j0,i_2,n-1))/(dT))^2 + 2*(PHI(j0,i_2,n-1) - PHI(j0,i_2,n-2))/dt -1))^(1/(gam-1));

                        % RHO(j-2, i) - (forward difference)
                        RHO_ij_1 = (1 - 0.5*(gam-1)*M0^2*(((PHI(j_1, i0, n-1) - PHI(j_2, i0, n2))/(dr))^2 + (1/RR(j_2,i0))^2*((PHI(j_2, i1, n2)-PHI(j_2,i_1,n2))/(2*dT))^2 + 2*(PHI(j_2,i0,n) - PHI(j_2,i0,n-1))/dt - 1))^(1/(gam-1));

                        % RHO_S(j,i)
                        PHI_Tij = (PHI(j0,i1, n-1) - PHI(j0,i_1, n-1))/(2*dT);
                        PHI_Rij = (PHI(j1,i0,n-1)-PHI(j_1,i0,n1))/(2*dr);
                        a2_ij = 2/(M0^2) - (gam-1) * (PHI_Rij^2 + ((1/RR(j0,i0))*PHI_Tij)^2 - 1 );
                        M_ij = (((1/RR(j0,i0))*PHI_Tij)^2 + PHI_Rij^2)/(a2_ij); % calculate local Mach number
                        e_ij = max([0.0, 1-(a2_ij/(M_ij^2))]); % artificial viscosity
                        RHO_Sij = e_ij*(PHI_Tij*(RHO_ij - RHO_i0j) + PHI_Rij*(RHO_ij - RHO_ij0));

                        % RHO_S(j+1,i)
                        PHI_Tij1 = (PHI(j1,i1, n-1) - PHI(j1,i_1, n-1))/(2*dT);
                        PHI_Rij1 = (PHI(j2,i0,n-1)-PHI(j0,i0,n-1))/(2*dr);
                        a2_ij1 = 2/(M0^2) - (gam-1) * (PHI_Rij1^2 + ((1/RR(j1,i0))*PHI_Tij1)^2 - 1 );
                        M_ij1 = (((1/RR(j1,i1))*PHI_Tij1)^2 + PHI_Rij1^2)/(a2_ij1);
                        e_ij1 = max([0.0, 1-(a2_ij1/(M_ij1^2))]);
                        RHO_Sij1 = e_ij1*(PHI_Tij1*(RHO_ij - RHO_i0j) + PHI_Rij1*(RHO_ij - RHO_ij0));

                        % RHO_S(j-1,i)
                        PHI_Tij0 = (PHI(j_1,i1, n1) - PHI(j_1,i_1, n1))/(2*dT);
                        PHI_Rij0 = (PHI(j0,i0,n-1)-PHI(j_2,i0,n2))/(2*dr);
                        a2_ij0 = 2/(M0^2) - (gam-1) * (PHI_Rij0^2 + ((1/RR(j_1,i0))*PHI_Tij0)^2 - 1 );
                        M_ij0 = (((1/RR(j_1,i0))*PHI_Tij0)^2 + PHI_Rij0^2)/(a2_ij0);
                        e_ij0 = max([0.0, 1-(a2_ij0/(M_ij0^2))]);
                        RHO_Sij0 = e_ij0*(PHI_Tij0*(RHO_ij0 - RHO_i0j0) + PHI_Rij0*(RHO_ij0 - RHO_ij_1)); % need j-2

                        % RHO_S(j,i+1), j = j0, i = i1
                        PHI_Ti1j = (PHI(j0,i2, n-1) - PHI(j0,i0, n-1))/(2*dT);
                        PHI_Ri1j = (PHI(j1,i1,n-1)-PHI(j_1,i1,n1))/(2*dr);
                        a2_i1j = 2/(M0^2) - (gam-1) * (PHI_Ri1j^2 + ((1/RR(j0,i1))*PHI_Ti1j)^2 - 1 );
                        M_i1j = (((1/RR(j0,i1))*PHI_Ti1j)^2 + PHI_Ri1j^2)/(a2_i1j);
                        e_i1j = max([0.0, 1-(a2_i1j/(M_i1j^2))]);
                        RHO_Si1j = e_i1j*(PHI_Ti1j*(RHO_i1j - RHO_ij) + PHI_Ri1j*(RHO_i1j - RHO_i1j0));

                        % RHO_S(j,i-1)
                        PHI_Ti0j = (PHI(j0,i0, n-1) - PHI(j0,i_2, n-1))/(2*dT);
                        PHI_Ri0j = (PHI(j1,i_1,n-1)-PHI(j_1,i_1,n1))/(2*dr);
                        a2_i0j = 2/(M0^2) - (gam-1) * (PHI_Ri0j^2 + ((1/RR(j0,i_1))*PHI_Ti0j)^2 - 1 );
                        M_i0j = (((1/RR(j0,i_1))*PHI_Ti0j)^2 + PHI_Ri0j^2)/(a2_i0j);
                        e_i0j = max([0.0, 1-(a2_i0j/(M_i0j^2))]);
                        RHO_Si0j = e_i0j*(PHI_Ti0j*(RHO_i0j - RHO_i_1j) + PHI_Ri0j*(RHO_i0j - RHO_i0j0)); % need RHO(i-2, j)

                        % Density averages
                        RHO_j0 = 0.5*(RHO_ij + RHO_ij0 - RHO_Sij - RHO_Sij0); % (i, j-0.5)
                        RHO_j1 = 0.5*(RHO_ij + RHO_ij1 - RHO_Sij - RHO_Sij1); % (i, j+0.5)
                        RHO_i0 = 0.5*(RHO_ij + RHO_i0j - RHO_Sij - RHO_Si0j); % (i-0.5, j)
                        RHO_i1 = 0.5*(RHO_ij + RHO_i1j - RHO_Sij - RHO_Si1j); % (i+0.5, j)

                        % Second Order Derivatives
                        PHI_TT = (1/RR(j0,i0)^2)*(1/dT)*(RHO_i1*(PHI(j0,i1,n-1)-PHI(j0,i0,n-1))/dT - RHO_i0*(PHI(j0,i0,n-1)-PHI(j0,i_1,n-1))/dT);
                        PHI_RR = (1/RR(j0,i0))*((0.5*(RR(j1,i0)+RR(j0,i0))*(PHI(j1,i0,n-1)-PHI(j0,i0,n-1))/dr)*RHO_j1 - RHO_j0*(0.5*(RR(j0,i0)+RR(j_1,i0))*(PHI(j0,i0,n-1)-PHI(j_1,i0,n-1))/dr))/dr;

                        % Apply Conditions for Time
                        PHI(j0,i0,n) = ((0.5*alpha/dt - 1/dt^2)*PHI(j0,i0,n-2) + 2*PHI(j0,i0,n-1)/dt^2 + PHI_RR + PHI_TT)/(1/dt^2 + 0.5*alpha/dt);
                    else
                        PHI(j,i,n) = PHI(j,i,n-1); % set the boundary condition
                    end
                end
            end
            difference = abs(PHI(:,:,n) - PHI(:,:,n-1));
            res(n-2) = max(difference(:));
        end

            %% Generate Plots

        % Uvals = size(PHI(:,

%         figure();
%         contourf(XX, YY, PHI(:,:,length(t_vals)));

        PHI_R = zeros(size(PHI));
        PHI_T = zeros(size(PHI));
        for n = 1:length(t_vals)
            for j = 2:length(r_vals)
               PHI_R(j,:,n) =  (PHI(j,:,n)-PHI(j-1,:,n))./(RR(j,:) - RR(j-1,:));
            end

            for i = 1:length(T_vals)
                i0 = i;
                if i == 1
                    i1 = i+1;
                    i_1 = length(T_vals);
                elseif i == length(T_vals)
                    i1 = 1;
                    i_1 = i-1;
                else
                    i1 = i+1;
                    i_1 = i-1;
                end 

                PHI_T(:,i0,n) =  (1/RR(:,i0))*0.5*(PHI(:,i1,n)-PHI(:,i_1,n))./(TT(:,i1) - TT(:,i_1));
            end
        end
        PHI_T1(mm,:) = PHI_T(1,:,end);
        PHI_RHO = (1 - 0.5.*(gam-1).*M0^2.*(PHI_R(:,:,end).^2 + (PHI_T(:,:,end)./RR(:,:)).^2 - 1)).^(1/(gam-1));
        PHI_RHO1(mm,:) = PHI_RHO(1,:,end);
    end
    res_time{tt} = [t_vals(1:end-2); res(1:end-2)];
end


figure();
contourf(XX, YY, 1-((PHI_R(:,:,end).^2 + PHI_T(:,:,end).^2)), 50); %./((RR.*cos(TT)).^2)
colorbar('eastoutside');

figure();
for mm = 1:length(M_vals)
    plot(T_vals, PHI_T1(mm,:));
    hold on;
end



% figure();
% quiver(XX,YY, PHI_)

% figure();
% contourf(XX, YY, (PHI_R(:,:,1).^2 + PHI_T(:,:,1).^2).^0.5);

%% Plot Density

% PHI_RHO = (1 - 0.5*(gam-1)*(PHI_R(:,:,end)^2 + (PHI_T(:,:,end)/RR(:,:))^2 - 1))^(1/(gam-1));
strtype = {'--', 'o'};
figure();
for mm = 1:length(M_vals)
    plot(T_vals, PHI_RHO1(mm,:), strtype{mm});
    hold on;
end

figure(); loglog(res_time{1}(1,:), res_time{1}(2,:));
hold on;
loglog(res_time{2}(1,:), res_time{2}(2,:), 'o');
hold on;
loglog(res_time{3}(1,:), res_time{3}(2,:), 'o');