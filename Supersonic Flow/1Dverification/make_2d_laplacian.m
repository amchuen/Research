function A = make_2d_laplacian(m)

n = m^2;

B = zeros(n,3);

B(:,1) = 2 * ones(n,1);

Xtmp = [2 : m]' * ones(1,m) + m * ones(m-1,1) * [0 : m-1];
I = reshape(Xtmp, n-m, 1);
B(I,2) = -ones(n-m,1);
B(:,3) = -ones(n,1);

A = spdiags(B,[0 1 m],n,n);

A = A + A';

