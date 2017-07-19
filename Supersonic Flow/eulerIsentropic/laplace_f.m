function laplacian = laplace_f(FF, dx, dy, BC, enforce)

%% laplace_f -> takes the laplacian of the given field variable FF

% This function returns the laplacian of the discrete FF. It calculates the
% second derivative in both directions of the field from the forward and 
% backward first-order, first derivatives at a given point IJ, and adds the
% values in both directions to return the numerical laplacian.

%%%% Inputs:
% * FF: field values
% * dx: spacing in the x-direction (column-wise)
% * dy: spacing in the y-direction (row-wise)
% * BC: struct input containing the neumann or dirichlet conditions 
% * enforce: enforces the derivative at the BC for central difference
%   calculation if it is Neumann

%%%% Other Notes:
% * Dirichlet conditions are enforced in the calculation through using
%   imaginary points one step-size outside the grid to determine the
%   derivatives at the boundary.


%% Solution

% Get F_xx
[~, fx_f, fx_b]= grad_f(FF, 2, dx, BC, enforce);
F_xx = (fx_f - fx_b)./dx;
 
% Get F_yy
[~, fy_f, fy_b]= grad_f(FF, 1, dy, BC, enforce);
F_yy = (fy_f - fy_b)./dy;
 
% Get Full Laplacian
laplacian = F_xx + F_yy;

end