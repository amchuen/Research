clc;
clear;
close all;

%% Define Grid Space and Control Points
dx = pi/8;
dy = dx;
x_start = 0;
x_end = 2*pi;
y_start = 0;
y_end = 2*pi;

res = 100; % resolution of the plots

xvals = linspace(x_start, x_end, (x_end - x_start)/dx + 1);
yvals = linspace(y_start, y_end, (y_end - y_start)/dy + 1);

[XX, YY] = meshgrid(xvals, yvals);
[y_ct, x_ct] = size(XX);

ZZ = sin(XX).*cos(YY);

% analytical derivatives
dfdx = @(x, y) cos(x).*cos(y);
dfdy = @(x,y) -sin(x).*sin(y);
dfdxy = @(x,y) -cos(x).*sin(y);

figure();
surf(XX, YY, ZZ);
% axis([0 dx 0 dy -3 3]); 
grid on;
hold on;
a_mat = zeros([size(XX), 4,4]);

%% Get Fx

fx = zeros(size(XX));
% a_mat(:,:,1) = ZZ;

for row = 1:y_ct
    [~, ~, a_x, b_x, c_x, d_x] = cubic_interp(XX(row,:), ZZ(row,:), res, dfdx(XX(row,1),YY(row,1)), dfdx(XX(row,end), YY(row,end)));
    
    fx(row,1:end) = [b_x, dfdx(XX(row,end), YY(row,end))];
%     a_mat(row, 1:end, 1, 2) = [b_x, dfdx(XX(row,end), YY(row,end))];
%     a_mat(row, 1:end-1, 1, 3) = c_x;
%     a_mat(row, 1:end-1, 1, 4) = d_x;
end

%% Get Fy

fy = zeros(size(XX));

for col = 1:x_ct
    [~, ~, a_y, b_y, c_y, d_y] = cubic_interp(YY(:,col)', ZZ(:,col)', res, dfdy(XX(1,col),YY(1,col)), dfdy(XX(end,col), YY(end, col)));
    
    fy(1:end, col) = [b_y'; dfdy(XX(end,col), YY(end, col))];
%     a_mat(1:end, col, 2, 1) = [b_y'; dfdy(XX(end,col), YY(end, col))];
%     a_mat(1:end-1, col, 3, 1) = c_y';
%     a_mat(1:end-1, col, 4, 1) = d_y';
end

%% Get Fxy

fxy = zeros(size(XX));

% Boundary conditions?
fxy(:,1) = dfdxy(XX(:,1),YY(:,1));
fxy(:,end) = dfdxy(XX(:,end),YY(:,end));
fxy(1,:) = dfdxy(XX(1,:),YY(1,:));
fxy(end,:) = dfdxy(XX(end,:),YY(end,:));


for row = 2:y_ct-1
    for col = 2:x_ct-1
        fxy(row, col) = (ZZ(row+1,col+1) + ZZ(row-1, col-1) - ZZ(row+1, col-1) - ZZ(row-1, col+1))/(4*dx*dy);        
    end
end

% a_mat(:,:, 2,2) = fxy;


%% Solve for the rest of the coefficients using Gauss-Siedel?

Abar = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        1, dx, dx^2, dx^3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        1, 0, 0, 0, dy, 0, 0, 0, dy^2, 0, 0, 0, dy^3, 0, 0, 0;...
        1, dx, dx^2, dx^3, dy, dx*dy, dx^2*dy, dx^3*dy, dy^2, dx*dy^2, dx^2*dy^2, dx^3*dy^2, dy^3, dx*dy^3, dx^2*dy^3, dx^3*dy^3;...
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 1, 2*dx, 3*dx^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 1, 0, 0, 0, dy, 0, 0, 0, dy^2, 0, 0, 0, dy^3, 0, 0;... 
        0, 1, 2*dx, 3*dx^2, 0, dy, 2*dx*dy, 3*dx^2*dy, 0, dy^2, 2*dx*dy^2, 3*dx^2*dy^2, 0, dy^3, 2*dx*dy^3, 3*dx^2*dy^3;... #fx(1,1)
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 1, dx, dx^2, dx^3, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 1, 0, 0, 0, 2*dy, 0, 0, 0, 3*dy^2, 0, 0, 0,;... 
        0, 0, 0, 0, 1, dx, dx^2, dx^3, 2*dy, 2*dx*dy, 2*dx^2*dy, 2*dx^3*dy, 3*dy^2, 3*dx*dy^2, dx^2*3*dy^2, dx^3*3*dy^2;... # fy(1,1)
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 0, 1, 2*dx, 3*dx^2, 0, 0, 0, 0, 0, 0, 0, 0;...
        0, 0, 0, 0, 0, 1, 0, 0, 0, 2*dy, 0, 0, 0, 3*dy^2, 0, 0, ;... 
        0, 0, 0, 0, 0, 1, 2*dx, 3*dx^2, 0, 2*dy, 2*2*dx*dy, 2*3*dx^2*dy, 0, 3*dx*dy^2, 2*dx*3*dy^2, 3*dx^2*3*dy^2;... # fxy(1,1)
        ];
    
        
for jj = 1:y_ct-1
    for ii = 1:x_ct-1
       fbar = [ZZ(jj, ii);
               ZZ(jj, ii+1);
               ZZ(jj+1, ii);
               ZZ(jj+1, ii+1);
               fx(jj, ii);
               fx(jj, ii+1);
               fx(jj+1, ii);
               fx(jj+1, ii+1);
               fy(jj, ii);
               fy(jj, ii+1);
               fy(jj+1, ii);
               fy(jj+1, ii+1);
               fxy(jj, ii);
               fxy(jj, ii+1);
               fxy(jj+1, ii);
               fxy(jj+1, ii+1);
               ];
        
        a_calc = Abar\fbar;
        a_mat(jj, ii, 1, 1) = a_calc(1);
        a_mat(jj, ii, 2, 1) = a_calc(2);
        a_mat(jj, ii, 3, 1) = a_calc(3);
        a_mat(jj, ii, 4, 1) = a_calc(4);
        a_mat(jj, ii, 1, 2) = a_calc(5);
        a_mat(jj, ii, 2, 2) = a_calc(6);
        a_mat(jj, ii, 3, 2) = a_calc(7);
        a_mat(jj, ii, 4, 2) = a_calc(8);
        a_mat(jj, ii, 1, 3) = a_calc(9);
        a_mat(jj, ii, 2, 3) = a_calc(10);
        a_mat(jj, ii, 3, 3) = a_calc(11);
        a_mat(jj, ii, 4, 3) = a_calc(12);
        a_mat(jj, ii, 1, 4) = a_calc(13);
        a_mat(jj, ii, 2, 4) = a_calc(14);
        a_mat(jj, ii, 3, 4) = a_calc(15);
        a_mat(jj, ii, 4, 4) = a_calc(16);
    end
end

%% Test Plot the first square?

% figure()
for i = 1:y_ct-1
    for ii = 1:x_ct-1
        xtest = linspace(XX(i,ii), XX(i,ii)+dx, 5);
        ytest = linspace(YY(i,ii), YY(i,ii)+dy, 5);

        [xx, yy] = meshgrid(xtest, ytest);
        zz = zeros(size(xx));

        for row = 1:4
            for col = 1:4
                zz = zz + a_mat(i,ii, col, row)*(((xx - XX(i,ii)).^(col-1))).*(((yy - YY(i,ii)).^(row-1)));
            end

        end

        plot3(xx,yy,zz, '.');
        hold on;
        
    end
end