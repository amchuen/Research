clear;
close all;
clc;

% Define the desired function and its derivative
f_x = @(x) tan(x) - 1;
dfdx = @(x) (sec(x)).^2;

% Make guess and desired tolerance
guess_init = 3.1/3.0;
tol = 1e-10; % Guess to 10 decimal places

% Solve using Newton-Raphson Method
[approx, iter, res] = newt_raph(f_x, dfdx, guess_init, tol);

fprintf('Approximate solution of pi:%0.10f\n', approx * 4.0);
fprintf('Actual solution of pi: %0.10f\n', pi);
fprintf('Num. of iterations needed:%i\n', iter);
fprintf('Residual: %0.20e\n', res(end));
fprintf('Comparison with actual Pi: %0.20f\n', approx *4.0 - pi);

figure()
plot(1:iter, res);
xlabel('Number of Iterations');
ylabel('Residual');

figure()
loglog((log(1:iter)), res);

%% sqrt of 2

f_x = @(x) (x.^2) - 2;
dfdx = @(x) 2.*x;

% Make guess and desired tolerance
guess_init = 1.4;
tol = 1e-10;

% Solve using Newton-Raphson Method
[approx, iter, res] = newt_raph(f_x, dfdx, guess_init, tol);

fprintf('Approximate solution of sqrt(2):%0.10f\n', approx);
fprintf('Actual solution of sqrt(2): %0.10f\n', sqrt(2));
fprintf('Num. of iterations needed:%i\n', iter);
fprintf('Residual: %0.20f\n', res(end));
fprintf('Comparison with actual Pi: %0.20f\n', approx - sqrt(2));

figure()
plot(1:iter, res);
xlabel('Number of Iterations');
ylabel('Residual');

figure()
loglog((log(1:iter)), res);