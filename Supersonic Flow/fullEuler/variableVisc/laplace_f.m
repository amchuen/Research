function laplacian = laplace_f(FF, dx, dy, BC, varargin)

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
enforce = 1;

if all(size(varargin{1}) == size(FF))
    % if true, then input must be additional info for polar calculations
    RR = varargin{1};
    is_polar = 1;
else
    is_polar = 0;
end

if is_polar
    % Get F_xx (theta direction)
    [fx, fx_f, fx_b]= grad_f(FF, 2, dx, BC, enforce);
    F_xx = ((fx_f - fx_b)./dx)./(RR.^2);

    % Get F_yy (radial direction)
    rr_n = [0.5.*(RR(2:end,:) + RR(1:(end-1),:)); 0.5.*((RR(end,:)+dy) + RR(end,:))];
    rr_s = 0.5.*([2.*RR(1,:); RR(2:end,:) + RR(1:(end-1),:)]);
    [~, fy_f, fy_b]= grad_f(FF, 1, dy, BC, enforce);
    F_yy = ((rr_n.*fy_f - rr_s.*fy_b)./dy)./(RR);
else
    % Get F_xx
    [fx, fx_f, fx_b]= grad_f(FF, 2, dx, BC, enforce);
    F_xx = (fx_f - fx_b)./dx;

    % Get F_yy
    [~, fy_f, fy_b]= grad_f(FF, 1, dy, BC, enforce);
    F_yy = (fy_f - fy_b)./dy;
end
 
% Get Full Laplacian
laplacian = F_xx + F_yy;

end