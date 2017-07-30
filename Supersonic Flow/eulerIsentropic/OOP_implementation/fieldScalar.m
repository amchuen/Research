classdef fieldScalar %< handle

%% ABOUT FIELD SCALAR
% generates field scalar values and methods to calculate derivs of field
% scalar
    
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

    properties
        % setup
        gr % grid
        fv % field value
        bc % boundary conditions
        ct % control variables

        % gradient derivs
        f_1f % gradient in row direction (1), forward (f)
        f_1b % gradient in row direction (1), backward (b)
        f_1 % central difference gradient

        f_2f
        f_2b
        f_2

        % laplacian
        f_11 % 2nd deriv in row direction
        f_22 % 2nd deriv in column direction
        laplace
       
    end
   
    methods
        function obj = fieldScalar(gr, fv, bc, ct)
           obj.gr = gr;
           obj.fv = fv;
           obj.bc = bc;
           obj.ct = ct;

        end

        function obj = Laplacian(obj)
            % Get f_11 and f_22
            if obj.ct.is_polar
                % deriv 1 is radial direction
                rr_n = [0.5.*(obj.gr.RR(2:end,:) + obj.gr.RR(1:(end-1),:)); 0.5.*((obj.gr.RR(end,:)+dr) + obj.gr.RR(end,:))];
                rr_s = 0.5.*([2.*obj.gr.RR(1,:); obj.gr.RR(2:end,:) + obj.gr.RR(1:(end-1),:)]);
                obj.f_11 = (rr_n.*obj.f_1f - rr_s.*obj.f_1b)./(obj.gr.d1.*obj.gr.RR);

                % deriv 2 is in angular direction
                obj.f_22 = (obj.f_2f - obj.f_2b)./(obj.gr.d2.*obj.gr.RR);

            else
                % deriv 1 is in y direction
                obj.f_11 = (obj.f_1f - obj.f_1b)./obj.gr.d1;

                % deriv 2 is in x direction
                obj.f_22 = (obj.f_2f - obj.f_2b)./obj.gr.d2;
            end

            % Get Full Laplacian
            obj.laplace= obj.f_11 + obj.f_22;

        end
        
        function obj = Gradient(obj) % maybe should create a vector instead?
            %% Solution
            if length(size(obj.fv)) == 3
                % means variable is being propagated using three-level
                % scheme
                FF = obj.fv(:,:,2);
            else
                % means variable is temporarily defined for current calc
                FF = obj.fv;
            end

            % Initialize Variables
            obj.f_1f = zeros(size(obj.fv));
            obj.f_1b = obj.f_1f;
            obj.f_1 = obj.f_1f;
            obj.f_2f = obj.f_1f;
            obj.f_2b = obj.f_1f;
            obj.f_2 = obj.f_1f;

            % Calculate forward and backward 1st-order differences
            obj.f_1f(1:(end-1),:) = (diff(FF,1,1)./obj.gr.d1);
            obj.f_1b(2:end,:) = diff(FF,1,1)./obj.gr.d1;
            if obj.ct.is_polar
                obj.f_2f(:,1:(end-1)) = (diff(FF,1,2)./(obj.gr.d2 .* obj.gr.RR(:,1:(end-1))));
                obj.f_2b(:,2:end) = diff(FF,1,2)./(obj.gr.d2 .* obj.gr.RR(:,2:end));
            else
                obj.f_2f(:,1:(end-1)) = (diff(FF,1,2)./obj.gr.d2);
                obj.f_2b(:,2:end) = diff(FF,1,2)./obj.gr.d2;
            end

            % control switches for enforcing neumann derivs if triggered
            ny = 0; sy = 0; ex = 0; wx = 0;

            % Loop through boundary conditions and apply as necessary
            fnames = fieldnames(obj.bc);
            for i = 1:length(fnames)
                if strcmp(fnames{i}, 'N') % use loop to check for 
                    obj.f_1f(end,:) = (obj.bc.N - FF(end,:))./obj.gr.d1;
                elseif strcmp(fnames{i}, 'Ny')
                    obj.f_1f(end,:) = obj.bc.Ny;
                    ny = 1;
                end

                if strcmp(fnames{i}, 'S')
                    % use reflective boundary condition if at wall
                    % i.e. scalars are equivalent accross boundaries,
                    obj.f_1b(1,:) = (FF(1,:) - BC.S)./obj.gr.d1;
                elseif strcmp(fnames{i}, 'Sy')
                    obj.f_1b(1,:) = BC.Sy; % assumes boundary body is between grid
                    sy = 1;
                end

                if obj.ct.is_polar
                    if strcmp(fnames{i}, 'E')
                        obj.f_2f(:,end) = (obj.bc.E - FF(:,end))./(obj.gr.d2.*obj.gr.RR(:,end));
                    elseif strcmp(fnames{i}, 'Ex')
                        f_exit = obj.bc.Ex .* (2.*obj.gr.d2.*obj.gr.RR(:,end)) + FF(:,end-1);
                        obj.f_2f(:,end) = (f_exit - FF(:,end))./(obj.gr.d.*obj.gr.RR(:,end));
                        ex = 1;
                    end

                    if strcmp(fnames{i}, 'W')
                        obj.f_2b(:,1) = (FF(:,1) - obj.bc.W)./(obj.gr.d2.*obj.gr.RR(:,end));
                    elseif strcmp(fnames{i}, 'Wx')
                        f_west = FF(:,2) - obj.bc.Wx.*2.*obj.gr.d2.*obj.gr.RR(:,end);
                        obj.f_2b(:,1) = (FF(:,1) - f_west)./(obj.gr.d2.*obj.gr.RR(:,end));
                        wx = 1;
                    end

                else
                    if strcmp(fnames{i}, 'E')
                        obj.f_2f(:,end) = (obj.bc.E - FF(:,end))./obj.gr.d2;
                    elseif strcmp(fnames{i}, 'Ex')
                        f_exit = obj.bc.Ex .* (2.*obj.gr.d2) + FF(:,end-1);
                        obj.f_2f(:,end) = (f_exit - FF(:,end))./obj.gr.d2;
                        ex = 1;
                    end

                    if strcmp(fnames{i}, 'W')
                        obj.f_2b(:,1) = (FF(:,1) - obj.bc.W)./obj.gr.d2;
                    elseif strcmp(fnames{i}, 'Wx')
                        f_west = FF(:,2) - obj.bc.Wx.*2.*obj.gr.d2;
                        obj.f_2b(:,1) = (FF(:,1) - f_west)./obj.gr.d2;
                        wx = 1;
                    end

                end

            end

            % Calculate Central Difference
            obj.f_1 = 0.5.*(obj.f_1f + obj.f_1b);
            obj.f_2 = 0.5.*(obj.f_2f + obj.f_2b);
            
            % enforce Neumann boundary conditions for central differences
            if ny
               obj.f_1(end,:) = obj.bc.Ny;
            end

            if sy
               obj.f_1(1,:) = obj.bc.Sy;
            end

            if ex
              obj.f_2(:,end) = obj.bc.Ex; 
            end

            if wx
               obj.f_2(:,1) = obj.bc.Wx;
            end


        end

        function obj = step_time(obj, varargin)
            % check for three-level
            if length(size(obj.fv))~=3
               error('This variable cannot be stepped forward in time!\n');
            end

            % calculate new "time-step"
            summation = zeros(size(ob.fv(:,:,3)));
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    summation = summation + varargin{i};
                end
            end
            obj.fv(:,:,3) = (obj.ct.eps_s.*obj.laplace - summation - (obj.ct.eps_t./(obj.ct.dt^2) - 0.5./obj.ct.dt).*obj.fv(:,:,1) + 2.*obj.ct.eps_t.*obj.fv(:,:,2)./(obj.ct.dt.^2))./(obj.ct.eps_t./(obj.ct.dt^2) + 0.5./obj.ct.dt);

            % update three-level
            obj.fv(:,:,1:2) = obj.fv(:,:,2:3);
            
            % update bounary conditions

       end
       
       % function setup_bc(obj)
       
   end
    
end