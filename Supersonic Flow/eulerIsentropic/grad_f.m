function [fx_c, fx_f, fx_b]= grad_f(FF, DIR, dx, BC, enforce_cent_BC)

%% grad_f -> takes the gradient of the given field variable FF in a desired direction

% This function returns the laplacian of the discrete FF. It calculates the
% second derivative in both directions of the field from the forward and 
% backward first-order, first derivatives at a given point IJ, and adds the
% values in both directions to return the numerical laplacian.

%%%% Inputs:
% * FF: field values
% * dx: spacing in the x-direction (can be single value, or sized equivalently to FF's size, specified thru DIR)
% * DIR: specifies axis/direction by which the derivative is taken
% * BC: struct input containing the neumann or dirichlet conditions 
% * enforce: enforces the derivative at the BC for central difference
%   calculation if it is Neumann

%%%% Other Notes:
% * Dirichlet conditions are enforced in the calculation through using
%   imaginary points one step-size outside the grid to determine the
%   derivatives at the boundary.

%% Solution

% Initialize Variables
fx_f = zeros(size(FF));
fx_b = fx_f;

% Calculate forward and backward 1st-order differences
if DIR == 1 % gradient is in the y-direction
    fx_f(1:(end-1),:) = (diff(FF,1,1)./dx);
    fx_b (2:end,:)= diff(FF,1,1)./dx;
    
    % control switches for enforcing neumann derivs if triggered
    ny = 0;
    sy = 0; 
    
    % Loop through boundary conditions and apply as necessary
    fnames = fieldnames(BC);
    for i = 1:length(fnames)
        if strcmp(fnames{i}, 'N') % use loop to check for 
            fx_f(end,:) = (BC.N - FF(end,:))./dx;
        elseif strcmp(fnames{i}, 'Ny')
            fx_f(end,:) = BC.Ny;
            ny = 1;
        end
        
        if strcmp(fnames{i}, 'S')
            % wall is usually assigned to this part of the field... want to
            % change step-size to represent wall in-between grid
            fx_b(1,:) = (FF(1,:) - BC.S)./(0.5*dx);
        elseif strcmp(fnames{i}, 'Sy')
            fx_b(1,:) = BC.Sy; % assumes boundary body is between grid
            sy = 1;
        end
        
    end
    
elseif DIR == 2 % gradient is in the x-direction ( or theta direction)
    
    if isscalar(dx) % convert to matrix to simplify calculations
        dx = dx.*ones(size(FF));
%         dx_b = dx.*ones(size(FF));
%     else
%         % expect input to be RR.*dT
%         dx_f = dx;
%         dx_b = dx;
        
    end
    
    fx_f(:,1:(end-1)) = (diff(FF,1,2)./dx(:,1:(end-1)));
    fx_b(:,2:end) = diff(FF,1,2)./dx(:,2:end);
    
    % control switches for enforcing neumann derivs if triggered
    ex = 0;
    wx = 0;
    
    % Loop through boundary conditions and apply as necessary
    fnames = fieldnames(BC);
    for i = 1:length(fnames)
        if strcmp(fnames{i}, 'E')
            fx_f(:,end) = (BC.E - FF(:,end))./ dx(:,end);
        elseif strcmp(fnames{i}, 'Ex')
            f_exit = BC.Ex .* (dx(:,end-1) + dx(:,end)) + FF(:,end-1);
%             fx_f(:,end) = BC.Ex;
            fx_f(:,end) = (f_exit - FF(:,end))./dx(:,end);
            ex = 1;
        end
        
        if strcmp(fnames{i}, 'W')
            fx_b(:,1) = (FF(:,1) - BC.W)./dx(:,1);
        elseif strcmp(fnames{i}, 'Wx')
            f_west = FF(:,2) - BC.Wx.*2.*dx(:,2);
            fx_b(:,1) = (FF(:,1) - f_west)./dx(:,1);
%             fx_b(:,1) = BC.Wx;
            wx = 1;
        end
        
    end 
end

% Calculate Central Difference
fx_c = 0.5.*(fx_f + fx_b);
if enforce_cent_BC % enforce BCs if necessary
   if DIR == 1
       if ny
           fx_c(end,:) = BC.Ny;
       end
       
       if sy
           fx_c(1,:) = BC.Sy;
       end
       
   elseif DIR == 2
       if ex
          fx_c(:,end) = BC.Ex; 
       end
       
       if wx
           fx_c(:,1) = BC.Wx;
       end
       
   end
end


end