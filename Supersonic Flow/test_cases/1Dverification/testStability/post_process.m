clc;
clear
close all;

%% Load Data

% load('ratio_1.mat', 'ratio_list');
load('ratio_1.mat', 'tol');

data = [];

for i = 1:11
    load(['ratio_' num2str(i) '.mat'], 'res');
    load(['ratio_' num2str(i) '.mat'], 'ratio');
    if (norm(res(end,:)) < tol) && (ratio ~= 0)
        data(end+1,:) = [ratio, length(res)];
    end
    clear res;
    clear ratio;
end

%% Plot Data

figure();
semilogy(data(:,1), data(:,2));
hold on;
semilogy(data(data(:,2)==min(data(:,2)),1), min(data(:,2)), 'o');
xlabel('$$\frac{\epsilon_x \Delta t^2}{\alpha \Delta x^2}$$','Interpreter','latex');
ylabel('No. iterations');
title('Convergence rate vs. Time-Viscosity to Spatial-Viscosity Ratio');
saveas(gcf, 'beta_ideal.pdf');