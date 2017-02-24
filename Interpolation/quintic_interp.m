function [xout, yout, cc, resid_k, resid_g, res1k, res1g] = quintic_interp(xin, yin, k1, kf, g1, gf, SOR, res)

% This function will generate a quintic interpolation of the inputted data
% points. Each spline generated takes the form of 

% S(x) = f_i + k_i(x - x_i) + 0.5*g_i (x - x_i)^2 + d_i (x - x_i)^3 + e_i
% (x - x_i)^4 + f_i (x - x_i)^5

% The function will iterate to generate k_i's and g_i values. The
% coefficients will then be passed into equations to generate d_i, e_i, and
% f_i values for all splines.

%% Setups
nvals = length(xin);
d_i = zeros(1, nvals);
e_i = zeros(1, nvals);
f_i = zeros(1, nvals);

k_i = zeros(1, nvals);
g_i = zeros(1, nvals);

% x spacing
C_x = diff(xin);
dy_i = diff(yin);

% use finite differences to get derivative guesses
k_i(1) = k1;
g_i(1) = g1;

for i = 2:length(xin)-1
    k_i(i) = (yin(i+1) - yin(i-1)) /(2*C_x(i));
    g_i(i) = (yin(i+1) - 2*yin(i) + yin(i-1))/(C_x(i)^2);
end

k_i(end) = kf;
g_i(end) = gf;


%% Iterate to get the set of k_i and g_i values

resid_k = 1.0;
resid_g = 1.0;
res1k = 0;
res1g = 0;

% initialize g-solvers
alf1 = zeros(size(xin));
bet1 = zeros(size(xin));
gam1 = zeros(size(xin));
del1 = zeros(size(xin));

% initialize k-solvers
alf2 = alf1;
bet2 = bet1;
gam2 = gam1;
del2 = del1;

cc = 0; % iteration counter

while (max(resid_k) >= 1e-5) || (max(resid_g) >= 1e-5)
    
    % equation 1 - solve g_i's
    alf1 = [0, -1./(4.*C_x(1:nvals-2)), 0];
    bet1 = [1, 0.75.*(1./C_x(1:nvals-2) + 1./C_x(2:nvals-1)), 1];
    gam1 = [0, -1./(4*C_x(2:nvals-1)), 0];
    del1 = [g1, 5*(dy_i(2:nvals-1)./(C_x(2:nvals-1).^3) - dy_i(1:nvals-2)./(C_x(1:nvals-2).^3))...
        + 2*k_i(1:nvals-2)./(C_x(1:nvals-2).^2)...
        + 3.*k_i(2:nvals-1).*(1./(C_x(1:nvals-2).^2) - 1./(C_x(2:nvals-1).^2))...
        - 2.*k_i(3:nvals)./(C_x(2:nvals-1).^2),...
        gf];
    %     
    g_star = thomas3(alf1, bet1, gam1, del1);

    % equation 2 - solve k_i's
    alf2 = [0, 7./(C_x(1:nvals-2).^3), 0];
    bet2 = [1, 8.*(1./(C_x(1:nvals-2).^3) + 1./(C_x(2:nvals-1).^3)), 1];
    gam2 = [0, 7./(C_x(2:nvals-1).^3), 0];
    del2 = [k1, 15.*(dy_i(2:nvals-1)./(C_x(2:nvals-1).^4) + dy_i(1:nvals-2)./(C_x(1:nvals-2).^4))...
        - g_i(1:nvals-2)./(C_x(1:nvals-2).^2)...
        + 1.5.*g_i(2:nvals-1).*(1./(C_x(1:nvals-2).^2) - 1./(C_x(2:nvals-1).^2))...
        + g_star(3:nvals)./(C_x(2:nvals-1).^2),...
        kf];

    k_star = thomas3(alf2, bet2, gam2, del2);

    % store old k_i vector
    k_old = k_i;
    g_old = g_i;

    % make new values with successive under relaxation
    k_i = SOR*(k_star - k_old) + k_old;
    g_i = SOR*(g_star - g_old) + g_old;

    % check for convergence
    if cc==0
       res1k = abs(k_i - k_old);
       res1g = abs(g_i - g_old);
    end
    resid_k = (abs(k_i - k_old));
    resid_g = (abs(g_i - g_old));
    cc = cc + 1;
        if cc > 10000
            fprintf('Iterations may not be converging.\n');
            break
        end
end

%% Generate the other set of constants
for i = 1:nvals-1
   d_i(i) = (1/(C_x(i)^2) ) * ( 10 * (yin(i+1)-yin(i))/C_x(i) - 6*k_i(i) - 4*k_i(i+1) - 3*C_x(i)*g_i(i)*0.5 + C_x(i) * g_i(i+1) * 0.5);
   e_i(i) = (1/(C_x(i)^3) ) * (7*k_i(i+1) + 8*k_i(i) - 15*(dy_i(i))/C_x(i) - C_x(i)*g_i(i+1) + 3*C_x(i)*g_i(i)*0.5);
   f_i(i) = (1/(C_x(i)^4) ) * (6*(dy_i(i))/C_x(i) - 3*k_i(i) - 3*k_i(i+1) + C_x(i)*g_i(i+1)*0.5 - C_x(i)*g_i(i)*0.5);

end

%% Generate interpolation points
xout = [];
yout = [];
for i = 1:nvals-1
    % generate x and y values
    xvals = linspace(xin(i), xin(i+1), res+1);
    dxvals = xvals - xin(i);
    yvals = yin(i) + k_i(i) * dxvals + 0.5*g_i(i)*dxvals.^2 + d_i(i)*dxvals.^3 + e_i(i)*dxvals.^4 + f_i(i)*dxvals.^5;

    % append to out values
    xout = [xout, xvals];
    yout = [yout, yvals];
end

end