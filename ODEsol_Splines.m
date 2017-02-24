function u_vals = ODEsol_Splines(x_vals, u_0, u_n)
nvals = length(x_vals);
dx = x_vals(2) - x_vals(1);

a_vals = [0, (1/6 - (1/(dx^2)))*ones(1,nvals-2), 0];
b_vals = [1, (4/6 + 2/(dx^2))*ones(1,nvals-2), 1];
c_vals = a_vals;
d_vals = [u_0, zeros(1, nvals-2), u_n];

u_vals = thomas3(a_vals, b_vals, c_vals, d_vals);

end