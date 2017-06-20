clc;
clear;
close all;

%% Load Airfoil Coordinates

foil = importdata('BACNLF.dat');

%% Generate Physical Domain

% size of computational space
m_pts = 99; % number of points in x dir
n_pts = 50; % number of points in y dir

% % Physical space boundary conditions
% xS_p = linspace(0, 2, m_pts);
% xN_p = linspace(-0.5, 2.5, m_pts);
% xW_p = linspace(0, -0.5, n_pts);
% xE_p = linspace(2, 2.5, n_pts);
% 
% mW = -0.9/(0 - -0.5);
% mE = 0.9/(2.5-2);
% 
% hN = 1.0;
% kN = 1.3;
% aN = -0.4/(1.5^2);
% 
% yS_p = zeros(size(xS_p));
% yN_p = aN.*(xN_p - hN).^2 + kN;
% yW_p = mW .* (xW_p - -0.5) + 0.9;
% yE_p = mE .* (xE_p - 2.0);

x_LE = find(foil.data(:,1) == min(foil.data(:,1)), 1);
rad = 5.0;
xS_p = linspace(1,0, 50);
yS_p = spline(foil.data(1:x_LE,1), foil.data(1:x_LE,2), xS_p);
yS_p = [yS_p, spline(foil.data((x_LE+1):end,1), foil.data((x_LE+1):end,2), fliplr(xS_p(1:end-1)))];
xS_p = [xS_p, fliplr(xS_p(1:end-1))];

xW_p = ones(1, n_pts) * xS_p(1);
xE_p = xW_p;
xN_p = rad.*cos(linspace(0.5*pi, 1.5*pi, m_pts)) + 1;

yN_p = rad.*sin(linspace(0.5*pi, 1.5*pi, m_pts));
yE_p = linspace(spline(foil.data(1:x_LE,1), foil.data(1:x_LE,2),xW_p(1)), rad, n_pts); 
yW_p = linspace(spline(foil.data((x_LE+1):end,1), foil.data((x_LE+1):end,2),xW_p(end)), -rad, n_pts);

% figure();
% plot(xS_p, yS_p);hold on;
% plot(xN_p, yN_p);
% plot(xW_p, yW_p);
% plot(xE_p, yE_p);

% Initialize physical grid
XX = zeros(n_pts, m_pts);
YY = XX;

% Apply Boundary Conditions
XX(:,1) = xW_p;
XX(:,end) = xE_p;
XX(1,:) = xS_p;
XX(end,:) = xN_p;

YY(:,end) = yW_p;
YY(:,1) = yE_p;
YY(1,:) = yS_p;
YY(end,:) = yN_p;

% re-initialize physical grid, sweep from left to right and bottom to top
for i = 2:(size(XX,2)-1)
%     XX(2:end-1,i) = xS_p(i);
%     YY(2:end-1,i) = yW_p(2:end-1);
    XX(:,i) = linspace(XX(1,i), XX(end,i), n_pts);
    YY(:,i) = linspace(YY(1,i), YY(end,i), n_pts);
end

xvals = linspace(1,3, 25);
y_u = linspace(0,rad, n_pts);
y_l = linspace(0,-rad, n_pts);

[XX_u, YY_u] = meshgrid(fliplr(xvals(2:end)), y_u);
[XX_l, YY_l] = meshgrid(xvals(2:end), y_l);

XX = horzcat(XX_u, XX);
YY = horzcat(YY_u, YY);

XX = horzcat(XX, XX_l);
YY = horzcat( YY, YY_l);

figure();
plot(XX,YY, 'b-', XX', YY', 'b-');

%% Solve Elliptic Equation

res = 1;
tol = 1e-5;

while res > tol
    XX_old = XX;
    YY_old = YY;
   
    for i = 2:(size(XX,2)-1)
        a_xx = [0;...
                ((0.5*(XX(2:end-1,i+1)-XX(2:end-1,i-1))).^2 + (0.5*(YY(2:end-1,i+1) - YY(2:end-1,i-1))).^2);...
                0];
        b_xx = [1;...
                -2.*((0.5*(XX(2:end-1,i+1)-XX(2:end-1,i-1))).^2 + (0.5*(YY(2:end-1,i+1) - YY(2:end-1,i-1))).^2 + (0.5*(XX(3:end,i)-XX(1:end-2,i))).^2 + (0.5*(YY(3:end,i) - YY(1:end-2,i))).^2);...
                1];
        c_xx = [0;...
                ((0.5*(XX(2:end-1,i+1)-XX(2:end-1,i-1))).^2 + (0.5*(YY(2:end-1,i+1) - YY(2:end-1,i-1))).^2);...
                0];
        d_xx = [XX(1, i);...
                -(XX(2:end-1, i-1) + XX(2:end-1, i+1)).*((0.5*(XX(3:end,i)-XX(1:end-2,i))).^2 + (0.5*(YY(3:end,i) - YY(1:end-2,i))).^2) + ...
                2.*(0.5.*(XX(2:end-1,i+1)-XX(2:end-1,i-1)).*0.5.*(XX(3:end,i)-XX(1:end-2,i)) + 0.5.*(YY(2:end-1,i+1) - YY(2:end-1,i-1)).*0.5.*(YY(3:end,i) - YY(1:end-2,i))).*0.25.*(XX(3:end,i+1) + XX(1:end-2, i-1) - XX(3:end, i-1) - XX(1:end-2, i+1));...
                XX(end, i)];
        
        XX(:,i) = thomas3(a_xx, b_xx, c_xx, d_xx);
        
%         a_yy = [0; ones(size(YY(2:end-1,i-1)));0];
%         b_yy = [1; -4*ones(length(YY(2:end-1, i)),1); 1];
%         c_yy = [0; ones(size(YY(2:end-1,i-1)));0];
        d_yy = [YY(1, i);...
                -(YY(2:end-1, i-1) + YY(2:end-1, i+1)).*((0.5*(XX(3:end,i)-XX(1:end-2,i))).^2 + (0.5*(YY(3:end,i) - YY(1:end-2,i))).^2) + ...
                2.*(0.5.*(XX(2:end-1,i+1)-XX(2:end-1,i-1)).*0.5.*(XX(3:end,i)-XX(1:end-2,i)) + 0.5.*(YY(2:end-1,i+1) - YY(2:end-1,i-1)).*0.5.*(YY(3:end,i) - YY(1:end-2,i))).*0.25.*(YY(3:end,i+1) + YY(1:end-2, i-1) - YY(3:end, i-1) - YY(1:end-2, i+1));...
                YY(end, i)];
        
        YY(:,i) = thomas3(a_xx, b_xx, c_xx, d_yy);
    end
    
    res_x = max(abs(XX(:) - XX_old(:)));
    res_y = max(abs(YY(:) - YY_old(:)));
    
    res(end+1) = max(res_x, res_y);
    
end

figure();
plot(XX,YY, 'b-', XX', YY', 'b-');
axis equal