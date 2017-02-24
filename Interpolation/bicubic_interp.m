function [xout, yout, zout] = bicubic_interp(XX, YY, ZZ, dfdx_1, dfdx_n, dfdy_1, dfdy_n, dfdxy_x1, dfdxy_xn, dfdxy_y1, dfdxy_yn, res)

%% SETUP

% Define coefficient matrix
a_mat = zeros([size(XX), 4,4]);

% Define grid dimensions
[y_ct, x_ct] = size(XX);
dx = XX(1, 2) - XX(1,1);
dy = YY(2,1) - YY(1,1);

%% GET FX

fx = zeros(size(XX));
% a_mat(:,:,1) = ZZ;

for row = 1:y_ct
    [~, ~, a_x, b_x, c_x, d_x] = cubic_interp(XX(row,:), ZZ(row,:), res, dfdx_1(row), dfdx_n(row));
    
    fx(row,1:end) = [b_x, dfdx_n(row)];
%     a_mat(row, 1:end, 1, 2) = [b_x, dfdx(XX(row,end), YY(row,end))];
%     a_mat(row, 1:end-1, 1, 3) = c_x;
%     a_mat(row, 1:end-1, 1, 4) = d_x;
end

%% GET FY

fy = zeros(size(XX));

for col = 1:x_ct
    [~, ~, a_y, b_y, c_y, d_y] = cubic_interp(YY(:,col)', ZZ(:,col)', res, dfdy_1(col), dfdy_n(col));
    
    fy(1:end, col) = [b_y'; dfdy_n(col)];
%     a_mat(1:end, col, 2, 1) = [b_y'; dfdy(XX(end,col), YY(end, col))];
%     a_mat(1:end-1, col, 3, 1) = c_y';
%     a_mat(1:end-1, col, 4, 1) = d_y';
end


%% GET FXY

fxy = zeros(size(XX));

% Boundary conditions?
fxy(:,1) = dfdxy_x1;
fxy(:,end) = dfdxy_xn;
fxy(1,:) = dfdxy_y1;
fxy(end,:) = dfdxy_yn;


for row = 2:y_ct-1
    for col = 2:x_ct-1
        fxy(row, col) = (ZZ(row+1,col+1) + ZZ(row-1, col-1) - ZZ(row+1, col-1) - ZZ(row-1, col+1))/(4*dx*dy);        
    end
end

%% SOLVE FOR REST OF THE COEFFICIENTS

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

%% GENERATE INTERPOLATED POINTS
xout = zeros(res*(y_ct-1), res*(x_ct-1));
yout = zeros(size(xout));
zout = zeros(size(xout));

for i = 1:y_ct-1
    for ii = 1:x_ct-1
        xtest = linspace(XX(i,ii), XX(i,ii)+dx, res);
        ytest = linspace(YY(i,ii), YY(i,ii)+dy, res);

        [xx, yy] = meshgrid(xtest, ytest);
        zz = zeros(size(xx));

        for row = 1:4
            for col = 1:4
                zz = zz + a_mat(i,ii, col, row)*(((xx - XX(i,ii)).^(col-1))).*(((yy - YY(i,ii)).^(row-1)));
            end

        end
        
        xout(((i-1)*res+1):(i*res), ((ii-1)*res+1):(ii*res)) = xx;
        yout(((i-1)*res+1):(i*res), ((ii-1)*res+1):(ii*res)) = yy;
        zout(((i-1)*res+1):(i*res), ((ii-1)*res+1):(ii*res)) = zz;
    end
end

end