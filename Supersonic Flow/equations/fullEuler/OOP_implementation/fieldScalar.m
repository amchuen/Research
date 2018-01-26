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
        function obj = fieldScalar(gr, bc, ct, varargin)
           obj.gr = gr;
           obj.bc = bc;
           obj.ct = ct;
           
           if ~isempty(varargin)
               obj.fv = varargin{1};
               
           else
               obj.fv = ones(size(obj.gr.d22));
                            
           end
           
           obj = obj.Gradient();
           obj = obj.Laplacian();           

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
            % Solution
%             if length(size(obj.fv)) == 3
%                 % means vector does not have time dimension
%                 % third dimension is not related to domain, but rather
%                 % storage for variables
%                 FF = obj.fv(:,:,:,2);
%             else
%                 % means variable does not have time aspect
%                 FF = obj.fv;
%             end
            FF = obj.fv;
            
            % Initialize Variables
            obj.f_1f = zeros(size(FF));
            obj.f_1b = obj.f_1f;
            obj.f_1 = obj.f_1f;
            obj.f_2f = obj.f_1f;
            obj.f_2b = obj.f_1f;
            obj.f_2 = obj.f_1f;

            % Calculate forward and backward 1st-order differences
            obj.f_1f(1:(end-1),:,:) = (diff(FF,1,1)./obj.gr.d1);
            obj.f_1b(2:end,:,:) = diff(FF,1,1)./obj.gr.d1;
            if obj.ct.is_polar
                obj.f_2f(:,1:(end-1),:) = (diff(FF,1,2)./(obj.gr.d2 .* obj.gr.RR(:,1:(end-1))));
                obj.f_2b(:,2:end,:) = diff(FF,1,2)./(obj.gr.d2 .* obj.gr.RR(:,2:end));
            else
                obj.f_2f(:,1:(end-1),:) = (diff(FF,1,2)./obj.gr.d2);
                obj.f_2b(:,2:end,:) = diff(FF,1,2)./obj.gr.d2;
            end

            % Loop through boundary conditions and apply
            fnames = fieldnames(obj.bc);
            for i = 1:length(fnames) % loop through directions
                % generate temporary array of values from
                % scalars or vectors
                if strcmp(obj.bc.(fnames{i}).physical, 'sym')
                    BC_vals = obj.gen_symBC(FF, fnames{i}, obj.bc.(fnames{i}).numerical);
                    
                elseif strcmp(obj.bc.(fnames{i}).physical, 'wall')
                    BC_vals = obj.gen_wallBC(FF, fnames{i}, obj.bc.(fnames{i}).numerical, obj.bc.(fnames{i}).val);
                    
                else
                    BC_vals = obj.gen_regBC(FF, fnames{i}, obj.bc.(fnames{i}).numerical, obj.bc.(fnames{i}).val);
                    
                end                       
                
                % assign the boundary conditions to the corresponding
                % gradient matrices
                switch fnames{i}
                    case 'W'
                        obj.f_2b(:,1,:) = BC_vals;
                    case 'E'
                        obj.f_2f(:,end,:) = BC_vals;
                    case 'S'
                        obj.f_1b(1,:,:) = BC_vals;
                    case 'N'
                        obj.f_1f(end,:,:) = BC_vals;
                end
            end

            % Calculate Central Difference
            obj.f_1 = 0.5.*(obj.f_1f + obj.f_1b);
            obj.f_2 = 0.5.*(obj.f_2f + obj.f_2b);


        end
       
       % function setup_bc(obj)
       
    end

    methods % three-level-scheme related operations
        
        function obj = update_fv(obj, updateVals)
        % this function updates the field values for the fieldVector
            obj.fv = updateVals{1};
            
            % re-calculate derivatives
            obj = obj.Gradient();
            obj = obj.Laplacian(); 
        end
        
        function obj = update_bc(obj, dir, updateVals)
           % need to somehow access boundary conditions easily at relevant
           % locations based on wall and etc.
           
            for ii = 1:length(updateVals)
                if ~isempty(updateVals{ii})
                    obj.bc.(dir).val{ii} = updateVals{ii};
                end
            end
            
        end
    end
    
    methods % helper functions for calculating the gradients
        
        function BC_vals = gen_symBC(obj, FF, dir, num_types)
            BC_vals = [];
            for ii = 1:size(obj.fv,3) % loop through vector "elements"
                % symmetry will likely appear right on the
                % boundary of the grid... maybe need to
                % account for displaced symmetry?? 
                switch dir
                    case 'W'
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, (FF(:,1,ii) - FF(:,2,ii))./obj.gr.d2);
                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(:,1,ii) + FF(:,2,ii))./obj.gr.d2);
                        end 

                    case 'E'
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, (FF(:,end-1,ii) - FF(:,end,ii))./obj.gr.d2);
                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (-FF(:,end-1,ii) - FF(:,end,ii))./obj.gr.d2);
                        end

                    case 'N'
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, (FF(end-1,:,ii) - FF(end,:,ii))./obj.gr.d1);
                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (-FF(end-1,:,ii) - FF(end,:,ii))./obj.gr.d1);
                        end

                    case 'S'
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, (FF(1,:,ii) - FF(2,:,ii))./obj.gr.d1);
                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(1,:,ii) + FF(2,:,ii))./obj.gr.d1);
                        end 
                end

            end
            
        end
        
        function BC_vals = gen_regBC(obj, FF, dir, num_types, vals)
            BC_vals = [];
            for ii = 1:size(obj.fv, 3)
                switch dir
                    case 'W'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(:,1)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            f_BC = FF(:,2,ii) - temp_BC.*2.*obj.gr.d2;
                            BC_vals = cat(3, BC_vals, (FF(:,1,ii) - f_BC)./obj.gr.d2);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(:,1,ii) - temp_BC)./obj.gr.d2);

                        end 

                    case 'E'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(:,end)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            f_BC = temp_BC.*2.*obj.gr.d2 + FF(:,end-1,ii);
                            BC_vals = cat(3, BC_vals, (f_BC - FF(:,end,ii))./obj.gr.d2);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (temp_BC - FF(:,end,ii))./obj.gr.d2);

                        end 
                        
                    case 'N'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(end,:)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            f_BC = temp_BC.*2.*obj.gr.d2 + FF(end-1,:,ii);
                            BC_vals = cat(3, BC_vals, (f_BC - FF(end,:,ii))./obj.gr.d1);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (temp_BC - FF(end,:,ii))./obj.gr.d1);

                        end 
                        
                    case 'S'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(1,:)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            f_BC = FF(2,:,ii) - temp_BC.*2.*obj.gr.d1;
                            BC_vals = cat(3, BC_vals, (FF(1,:,ii) - f_BC)./obj.gr.d1);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(1,:,ii) - temp_BC)./obj.gr.d1);

                        end 
                        
                end
           end
        end
        
        function BC_vals = gen_wallBC(obj, FF, dir, num_types, vals)
            BC_vals = [];
            for ii = 1:size(obj.fv, 3)
                switch dir
                    case 'W'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(:,1)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, temp_BC);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(:,1,ii) - temp_BC)./(0.5.*obj.gr.d2));

                        end 

                    case 'E'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(:,end)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, temp_BC);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (temp_BC - FF(:,end,ii))./(0.5.*obj.gr.d2));

                        end 
                        
                    case 'N'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(end,:)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, temp_BC);

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (temp_BC - FF(end,:,ii))./(0.5.*obj.gr.d1));

                        end 
                        
                    case 'S'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(obj.gr.d11(1,:)));
                        else
                            temp_BC = vals{ii};
                        end
                        
                        % calculate boundary deriv
                        if strcmp(num_types{ii}, 'N')
                            BC_vals = cat(3, BC_vals, temp_BC); % need to switch order, or else derivs come out wrong!

                        elseif strcmp(num_types{ii}, 'D')
                            BC_vals = cat(3, BC_vals, (FF(1,:,ii) - temp_BC)./(0.5.*obj.gr.d1));

                        end 
                        
                end
            end
        end
    end
    
end