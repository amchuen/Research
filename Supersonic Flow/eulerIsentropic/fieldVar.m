classdef fieldVar
   properties
       GR % grid
   end
   
   methods
       function laplacian = laplace_f(FF, dx, dy, BC, enforce)

            % Get F_xx
            [~, fx_f, fx_b]= grad_f(FF, 2, dx, BC, enforce);
            F_xx = (fx_f - fx_b)./dx;

            % Get F_yy
            [~, fy_f, fy_b]= grad_f(FF, 1, dy, BC, enforce);
            F_yy = (fy_f - fy_b)./dy;

            % Get Full Laplacian
            laplacian = F_xx + F_yy;

        end
       
   end
    
end