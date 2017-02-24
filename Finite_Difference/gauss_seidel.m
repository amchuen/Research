function [u_guess, res, cc1] = gauss_seidel(Amat, bvals, W_load, dx, u_xx1, u_xxn, w, u_in, tol)
nvals = length(bvals);
u_guess = u_in;

res = max(abs((Amat * u_in')' - bvals));
cc1 = 1;

while res > tol
    
    for i = 1:nvals
       sigma = 0;
       
       for j = 1:nvals
           if j ~= i
               sigma = sigma + Amat(i,j)*u_guess(j);
           end
       end

       u_guess(i) = w*(bvals(i) - sigma)/Amat(i,i) + (1-w)*u_guess(i);
    end
    
    bvals = [0, -W_load(2) - (1/(dx^4))*(u_xx1 + u_guess(2) - 2*u_guess(3) + u_guess(4)), -W_load(3:end-2) - (1/(dx^4))*(u_guess(1:end-4) - 2*u_guess(2:end-3) + 2*u_guess(3:end-2) - 2*u_guess(4:end-1) + u_guess(5:end)), -W_load(end-1) - (1/(dx^4))*(u_xxn + u_guess(end-1) - 2*u_guess(end-2) + u_guess(end-3)), 0];
    
    res = max(abs((Amat * u_guess')' - bvals));
    cc1 = cc1 + 1;
end