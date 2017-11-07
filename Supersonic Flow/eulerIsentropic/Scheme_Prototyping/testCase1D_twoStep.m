clc;
close all;
clear;

%% Declare simulation parameters
dx = 2/400;
x_vals = (-1+dx):dx:(1-dx);
eps_arr = 1e-9;%[0.1, 0.01, 0.001];
dt = dx/10;
CFL_max = 0.95;
tol = 1e-5;
uu_1 = 1;
uu_n = -1;

u_plots = zeros(size(eps_arr));
res_plots = u_plots;
labels = {};

%% Initial Values and Boundary conditions
for i = 1:length(eps_arr)
    eps = eps_arr(i);
    labels{i} = ['eps = ' num2str(eps)];
    u_vals = repmat(-1.*x_vals, 3,1);
%     u_vals(:,(x_vals < 0.5)) = 1;
%     u_vals(:,(x_vals >= 0.5)) = -1;
    u_star = u_vals;


    %% Iterate
    res = [0];
    res_time = [0];
    while (res(end) > tol) ||  (length(res) < 300)
        u_vals(1:2,:) = u_vals(2:3,:);
        u_star(1:2,:) = u_star(2:3,:);

        CFL_i = dt.*max(abs(u_vals(2,:)))./(dx);

        if CFL_i > CFL_max
           fprintf('CFL condition not met!\n');
           fprintf('Decreasing time steps!\n');
           dt = dt * CFL_max / CFL_i;
           fprintf('New time step:%0.5f\n', dt);
        end

        % laplacian
        u_f = [u_star(2,2:end), uu_n];
        u_b = [uu_1, u_star(2,1:end-1)];
        laplace = (u_f + u_b)./dx^2;

        % spacial acceleration
        u2_xf = 0.5.*[diff(u_star(2,:).^2, 1,2), uu_n.^2 - u_star(2,end).^2]./dx;
        u2_xb = 0.5.*[u_star(2,1).^2 - uu_1.^2, diff(u_star(2,:).^2, 1,2)]./dx;
        u2_x = 0.5.*(u2_xf + u2_xb);
        
        % Second Step
        u_star(3,:) = ( eps.*laplace + u_star(1,:).*(1./(2*dt) - eps/(dx^2)) - u2_x )./(1./(2*dt) + eps./(dx^2));
        f_f = [diff(u_star([1 3],:), 1,2), uu_n - u_star([1 3], end)]./dx;
        f_b = [u_star([1 3], 1) - uu_1, diff(u_star([1 3],:),1,2)]./dx;
        f_vals = 0.5.*(f_f + f_b);
        
        u_vals(3,:) = u_vals(1,:) + 2.*dt.*(f_vals(1,:) + f_vals(2,:))./2;
        
        if all(isnan(u_vals(3,:)))
            fprintf('Solution failed!\n');
            fprintf('Iter:\t%0.5f\n\n', res(end));
            break;
        end

        if (size(res,1) == 1) && all(res(end,:) == 0)
            res(end) = max(abs((u_vals(3,:) - u_vals(2,:))));
            res_time(end) = dt;
        else
            res(end+1) = max(abs(u_vals(3,:) - u_vals(2,:)));
            res_time(end+1) = dt + res_time(end);
        end

    end

    %% Post Process
    figure(1);
    u_plots(i) = plot([-1, x_vals, 1], [1, u_vals(3,:), -1]);%, 'o-', 'MarkerIndices', 1:10:length([-1, x_vals, 1]));
    ylim([-1.5, 1.5]);
    title(['U(x,t), dx = ' num2str(dx)]);
    xlabel('x');
    ylabel('U');
    legend(u_plots(1:i), labels(1:i));
    hold on;

    figure(2);
    res_plots(i) = semilogy(res_time, res);
    title('Residual');
    xlabel('Pseudo-Time');
    ylabel('Residual (Error)');
    legend(res_plots(1:i), labels(1:i));
    hold on;
    
    pause(1)
end