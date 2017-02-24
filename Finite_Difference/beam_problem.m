clear;
close all;
clc;

%% Setup loads
nvals = 101;
xvals = linspace(0, 1, nvals);
dx = xvals(2) - xvals(1);
W_load = 1.*[0, ones(1,nvals-2), 0];
u_xx1 = 0;
u_xxn = u_xx1;

%% Analytical Solution
u_analytical = @ (x) -(x.^4 ./24) + x.^3.*xvals(end) ./12 - xvals(end)^3 .* x./24; 

%% Straight Shot Through PentaDiagonal Solver

E_pent = [(1/(dx^4)).*ones(1, nvals-3), 0];
A_pent = [-4/(dx^4), (-4/(dx^4)).*ones(1, nvals-3), 0];
D_pent = [1, 5/(dx^4), (6/(dx^4)).*ones(1, nvals - 4), 5/(dx^4), 1];
C_pent = [0, (-4/(dx^4)).*ones(1, nvals-3), -4/(dx^4), 0];
F_pent = [0, (1/(dx^4)).*ones(1, nvals-3), 0, 0];

B_pent = [0, -W_load(2:end-1), 0];

uvals_pent = pentadiagonal(E_pent, A_pent, D_pent, C_pent, F_pent, B_pent);

A_check = diag(E_pent, -2) + diag(A_pent, -1) + diag(D_pent) + diag(C_pent(1:end-1), 1) + diag(F_pent(1:end-2), 2);
test = A_check *uvals_pent';

figure();
plot(xvals, u_analytical(xvals), 'o');
hold on;
plot(xvals, uvals_pent);

%% Use the Pentadiagonal Lag Solution (Block Jacobi?)
sor_pent = 1.0;

% Run through first guess
u_pent = randi(4, size(xvals));
a_pent = [0, -2/(dx^4) .* ones(1, length(xvals)-2), 0];
b_pent = [1, 4/(dx^4) .* ones(1, length(xvals) -2), 1];
c_pent = [0, -2./(dx^4) .* ones(1, length(xvals)-2), 0];
% d_pent = [0, -W_load(2) - (1/(dx^4)).*(-2*u_pent(1) + u_pent(2) - 2*u_pent(3) + u_pent(4)),...
%     -W_load(3:end-2) - (1/(dx^4)).*(u_pent(1:end-4) - 2.*u_pent(2:end-3) + 2*u_pent(3:end-2) - 2*u_pent(4:end-1) + u_pent(5:end)),...
%     -W_load(end-1) - (1/(dx^4))*(-2*u_pent(end) + u_pent(end-1) - 2*u_pent(end-2) + u_pent(end-3)), 0];   
d_pent = zeros(size(xvals));
for i = 1:length(d_pent)
    if (i == 1) || (i == length(d_pent))
        d_pent(i) = 0;
    elseif (i == 2) 
        d_pent(i) = -W_load(i) - (1/(dx^4)).*(-2*u_pent(i-1) + u_pent(i) - 2*u_pent(i+1) + u_pent(i+2));
    elseif (i == (length(d_pent)-1))
        d_pent(i) = -W_load(i) - (1/(dx^4)).*(-2*u_pent(i+1) + u_pent(i) - 2*u_pent(i-1) + u_pent(i-2));
    else
        d_pent(i) = -W_load(i) - (1/(dx^4)).*(u_pent(i-2) - 2*u_pent(i-1) + 2*u_pent(i) - 2*u_pent(i+1) + u_pent(i+2));
    end
end


u_star = thomas3(a_pent, b_pent, c_pent, d_pent);
res = max(abs(u_star - u_pent));
cc = 0; % counter

while res > 1e-5
    % calculate the new guess
    u_pent = u_pent + sor_pent*(u_star - u_pent);
    
    % update the load
%     a_pent = [0, -2/(dx^4), -4./(dx^4) .* ones(1, length(xvals)-3), 0];
%     b_pent = [1, 5./(dx^4), 6/(dx^4) .* ones(1, length(xvals) -4), 5/(dx^4), 1];
%     c_pent = [0, -4./(dx^4) .* ones(1, length(xvals)-3), -2/(dx^4), 0];
    for i = 1:length(d_pent)
        if (i == 1) || (i == length(d_pent))
            d_pent(i) = 0;
        elseif (i == 2) 
            d_pent(i) = -W_load(i) - (1/(dx^4)).*(-2*u_pent(i-1) + u_pent(i) - 2*u_pent(i+1) + u_pent(i+2));
        elseif (i == (length(d_pent)-1))
            d_pent(i) = -W_load(i) - (1/(dx^4)).*(-2*u_pent(i+1) + u_pent(i) - 2*u_pent(i-1) + u_pent(i-2));
        else
            d_pent(i) = -W_load(i) - (1/(dx^4)).*(u_pent(i-2) - 2*u_pent(i-1) + 2*u_pent(i) - 2*u_pent(i+1) + u_pent(i+2));
        end
    end

    % calculate the new intermediate result
    u_star = thomas3(a_pent, b_pent, c_pent, d_pent);
    
    % generate residual
    res = max(abs(u_star - u_pent));
    cc = cc + 1;
end

hold on;
plot(xvals, u_pent);
hold on;

fprintf('Block Jacobi\n');
fprintf('Number of Iterations: %i\n', cc);
fprintf('Residual: %0.5e\n', res)
fprintf('Max Error with Analytical: %0.5e\n', max(abs(u_pent - uvals_pent)));

%% Use 2 explicit points
sor_exp = 1;
u_guess = randi(4, size(xvals));
u_guess(1) = 0;
u_guess(2) = u_analytical(xvals(2));
u_guess(end) = 0;
u_guess(end-1) = u_analytical(xvals(end-1));

a_exp = [0, 0, -2/(dx^4) .* ones(1, length(xvals)-4), 0, 0];
c_exp = [0, 0, -2./(dx^4) .* ones(1,length(xvals)-4), 0, 0];
b_exp = [1, 1,  4/(dx^4) .* ones(1, length(xvals)-4), 1, 1];

d_exp = zeros(size(xvals));
for i = 1:length(d_exp)
    if (i == 1) || (i == length(d_exp))
        d_exp(i) = 0;
    elseif (i == 2) || (i == (length(d_exp)-1))
        d_exp(i) = u_analytical(xvals(i));
    else
        d_exp(i) = -W_load(i) - (1/(dx^4)).*(u_guess(i-2) - 2*u_guess(i-1) + 2*u_guess(i) - 2*u_guess(i+1) + u_guess(i+2));
    end
% [0, u_analytical(xvals(2)),...
%         -W_load(3:end-2) - (1/(dx^4)).*(u_guess(1:end-4) - 2*u_guess(2:end-3) + 2*u_guess(3:end-2) - 2*u_guess(4:end-1) + u_guess(5:end)),...
%         u_analytical(xvals(end-1)), 0];
    
end

u_star = thomas3(a_exp, b_exp, c_exp, d_exp);
res1 = max(abs(u_star - u_guess));
cc1 = 0; % counter

while res1 > 1e-5
    % calculate the new guess
    u_guess = u_guess + sor_exp*(u_star - u_guess);
    
    % update the load
%     a_pent = [0, -2/(dx^4), -4./(dx^4) .* ones(1, length(xvals)-3), 0];
%     b_pent = [1, 5./(dx^4), 6/(dx^4) .* ones(1, length(xvals) -4), 5/(dx^4), 1];
%     c_pent = [0, -4./(dx^4) .* ones(1, length(xvals)-3), -2/(dx^4), 0];
%     d_exp = [0, u_analytical(xvals(2)),...
%         -W_load(3:end-2) - (1/(dx^4)).*(u_guess(1:end-4) - 2*u_guess(2:end-3) + 2*u_guess(3:end-2) - 2*u_guess(4:end-1) + u_guess(5:end)),...
%         u_analytical(xvals(end-1)), 0];
    for i = 1:length(d_exp)
        if (i == 1) || (i == length(d_exp))
            d_exp(i) = 0;
        elseif (i == 2) || (i == (length(d_exp)-1))
            d_exp(i) = u_analytical(xvals(i));
        else
            d_exp(i) = -W_load(i) - (1/(dx^4)).*(u_guess(i-2) - 2*u_guess(i-1) + 2*u_guess(i) - 2*u_guess(i+1) + u_guess(i+2));
        end
    end
    % calculate the new intermediate result
    u_star = thomas3(a_exp, b_exp, c_exp, d_exp);
    
    % generate residual
    res1 = max(abs(u_star - u_guess));
    cc1 = cc1 + 1;
end

fprintf('Explicit BC''s\n');
fprintf('Number of Iterations: %i\n', cc1);
fprintf('Residual: %0.5e\n', res1)
fprintf('Max Error with Analytical: %0.5e\n', max(abs(u_guess - uvals_pent)));

plot(xvals, u_guess);
legend('Analytical', 'Pent Solver', 'Lag Iteration', 'Explicit');
title('Validation of Various Numerical Methods for Solving Beams');
xlabel('Length (in)');
ylabel('Vertical Deflection (in)');
saveas(gcf, 'beam_validation.png');

% w = 0.35;
% tol = 1e-5;
% Amat = diag(a_pent(2:end), -1) + diag(b_pent, 0) + diag(c_pent(1:end-1), 1);
% u_in = ones(size(b_pent));
% bvals = generate_bvals(W_load, u_xx1, u_xxn, u_in, dx);
% 
% [u_SOR, resSOR, cc1] = gauss_seidel(Amat, bvals, W_load, dx, u_xx1, u_xxn, w, u_in, tol);
% 
% fprintf('Point Gauss Seidel\n');
% fprintf('Number of Iterations: %i\n', cc1);
% fprintf('Residual: %0.5e\n', res)