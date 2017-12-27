function y = rk4(tVec, func, y0)
    % RK4
    n = length(tVec);
    m = length(y0);
    y = zeros(m, n);
    dt = tVec(2) - tVec(1);
    y(:,1) = y0;

    for i = 1 : n - 1
        k1 = dt * func(y(:,i), tVec(i));
        k2 = dt * func(y(:,i) + k1 / 2, tVec(i) + dt / 2);
        k3 = dt * func(y(:,i) + k2 / 2, tVec(i) + dt / 2);
        k4 = dt * func(y(:,i) + k3, tVec(i) + dt);
        y(:,i+1) = y(:,i) + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end
end