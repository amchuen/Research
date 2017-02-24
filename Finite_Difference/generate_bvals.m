function dvals = generate_bvals(W_load, u_xx1, u_xxn, u_guess, dx)

dvals = [0, -W_load(2) - (1/(dx^4))*(u_xx1 + u_guess(2) - 2*u_guess(3) + u_guess(4)), -W_load(3:end-2) - (1/(dx^4))*(u_guess(1:end-4) - 2*u_guess(2:end-3) + 2*u_guess(3:end-2) - 2*u_guess(4:end-1) + u_guess(5:end)), -W_load(end-1) - (1/(dx^4))*(u_xxn + u_guess(end-1) - 2*u_guess(end-2) + u_guess(end-3)), 0];

end