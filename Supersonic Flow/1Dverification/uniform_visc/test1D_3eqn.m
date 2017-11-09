clc;
close all;
clear;

%% Declare simulation parameters
dx = 2/400;
x_vals = -1:dx:1;
eps_arr = [0.1, 0.01, 0.001];
dt = dx*0.01;
CFL_max = 0.95;
tol = 1e-5;

%% Initialize Values
i = 1;
eps = eps_arr(i);
labels{i} = ['eps = ' num2str(eps)];
fVals = ones(3,length(x_vals),4); %1:time   2:spatial   3:variables
fxVals = zeros(size(fVals(1,:,:)));

%% Boundary conditions
gam = 1.4;
M0 = 1.2;
P_l = 1./(gam .* M0^2);
rhoU_l = 1;
rho_l = 1;
rhoE_l = 1/(gam-1) * P_l + 0.5*rhoU_l^2;
fVals(1,1,:) = [rho_l;rhoU_l;rhoE_l; P_l];
fVals(1,end,:) = reshape([1,1,1,1],1,1,4).*fVals(1,1,:).*0.5;
fVals(1,end,end) = 1./(gam * (0.5*M0)^2);
for i = 1:size(fVals,3)
    fVals(1,:,i) = linspace(fVals(1,1,i), fVals(1,end,i), length(x_vals));
end
fVals(1,:,end) = (gam - 1) .* (fVals(1,:,3) - 0.5.*(fVals(1,:,2).^2 ./ fVals(1,:,1)));
fVals(2:3,:,:) = repmat(fVals(1,:,:),2,1,1);

%% Iterate
res = [0];
res_time = [0];
while (res(end) > tol) ||  (length(res) < 300)
    % Update fVals
    fVals(1:2,:,:) = fVals(2:3,:,:);
    fVals(2,:,end) = (gam - 1) .* (fVals(end,:,3) - 0.5.*(fVals(end,:,2).^2 ./ fVals(end,:,1))); % Pressure update
    fxVals = fVals(2,:,1:3);
    fxVals(:,:,3) = fxVals(:,:,3) + fVals(2,:,end); % add pressure to get enthalpy
    fxVals = fxVals .* repmat(fVals(2,:,2)./fVals(2,:,1),1,1,3);
    fxVals(:,:,2) = fxVals(:,:,2) + fVals(2,:,end); % add pressure to momentum
    
    CFL_i = dt.*max(abs(fVals(end,:,2)./fVals(end,:,1)))./dx;

    if CFL_i > CFL_max
       fprintf('CFL condition not met!\n');
       fprintf('Decreasing time steps!\n');
       dt = dt * CFL_max / CFL_i;
       fprintf('New time step:%0.5f\n', dt);
    end

    % laplacian
    laplace = (fVals(2,1:end-2,1:3) + fVals(2,3:end,1:3))./dx^2; % - 2.*fVals(2,2:end-1,1:3)

    % calculate advection/flux balance
    advec = (fxVals(:,3:end,:) - fxVals(:,1:end-2,:))./(2*dx);

    % Step Forward
    fVals(3,2:end-1,1:3) = (-advec + eps.*laplace + (0.5/dt - eps./(dx^2)).*fVals(1,2:end-1,1:3))./(0.5/dt + eps./(dx^2));
    
    if all(isnan(fVals(3,:)))
        fprintf('Solution failed!\n');
        fprintf('Iter:\t%0.5f\n\n', res(end));
        break;
    end

    if (size(res,1) == 1) && all(res(end,:) == 0)
        res(end) = max(abs((fVals(3,:) - fVals(2,:))));
        res_time(end) = dt;
    else
        res(end+1) = max(abs(fVals(3,:) - fVals(2,:)));
        res_time(end+1) = dt + res_time(end);
    end

end

%% Post Process
figure(1);
plot(x_vals, fVals(end,:,2));%, 'o-', 'MarkerIndices', 1:10:length([-1, x_vals, 1]));
title(['U(x,t), dx = ' num2str(dx)]);
xlabel('x');
ylabel('U');
% legend(u_plots(1:i), labels(1:i));

figure();
plot(x_vals, fVals(end,:,4));%, 'o-', 'MarkerIndices', 1:10:length([-1, x_vals, 1]));
title(['p(x,t), dx = ' num2str(dx)]);
xlabel('x');
ylabel('U');
% legend(u_plots(1:i), labels(1:i));
hold on;

figure();
semilogy(res_time, res);
title('Residual');
xlabel('Pseudo-Time');
ylabel('Residual (Error)');
hold on;