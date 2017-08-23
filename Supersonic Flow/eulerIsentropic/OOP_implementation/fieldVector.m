classdef fieldVector
    
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
    
    methods % initialization and calculus-based operations
        function obj = fieldVector(gr, bc, ct, time_axis, init_fv)
           obj.gr = gr;
           obj.bc = bc;
           obj.ct = ct;
           
           fv = [];
           for i = 1:length(init_fv)
              fv = cat(3, fv, init_fv{i});
           end
           
           if time_axis
               obj.fv = repmat(fv, 1, 1, 1, 3);...(-(1 + (r_cyl^2)./(obj.GR.d11.^2)).*sin(obj.GR.d22))), 1, 1, 1, 3); % rho * v
           else
               obj.fv = fv;                            
           end
           
           obj = obj.Gradient();
           obj = obj.Laplacian();           

        end
        
        function obj = Laplacian(obj)
            % Get f_11 and f_22
            if obj.ct.is_polar
                % deriv 1 is radial direction
                rr_n = [0.5.*(obj.gr.d11(2:end,:) + obj.gr.d11(1:(end-1),:)); 0.5.*((obj.gr.d11(end,:)+dr) + obj.gr.d11(end,:))];
                rr_s = 0.5.*([2.*obj.gr.d11(1,:); obj.gr.d11(2:end,:) + obj.gr.d11(1:(end-1),:)]);
                obj.f_11 = (rr_n.*obj.f_1f - rr_s.*obj.f_1b)./(obj.gr.d1.*obj.gr.d11);

                % deriv 2 is in angular direction
                obj.f_22 = (obj.f_2f - obj.f_2b)./(obj.gr.d2.*obj.gr.d11);

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
            if length(size(obj.fv)) == 4
                % means vector does not have time dimension
                % third dimension is not related to domain, but rather
                % storage for variables
                FF = obj.fv(:,:,:,2);
            else
                % means variable does not have time aspect
                FF = obj.fv;
            end
            
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
                BC_vals = obj.getBCvals(FF, fnames{i});
                
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
        
    end
    
    methods % three-level-scheme related operations
        
        function obj = update_fv(obj, updateVals)
        % this function updates the field values for the fieldVector
            if ((nargin == 1) && (length(size(obj.fv)) == 4))
               % vector has time dimension in it, and only needs to step
               % forward in time
               obj.fv(:,:,:,1:end-1) = obj.fv(:,:,:,2:end);
            elseif ((nargin == 2))
                % updated values must be passed into the new vector
                for i = 1:length(updateVals)
                   obj.fv(:,:,i) = updateVals{i}; 
                end 
            
            else
                error('Three-level scheme implementation error!\n Fix implementation or create a new update case.\n');
            end            
        end
        
        function obj = update_bc(obj, dir, bcIdx, updateVals)
           % need to somehow access boundary conditions easily at relevant
           % locations based on wall and etc.
           
            for ii = 1:length(updateVals)
                if ~isempty(updateVals{ii})
                    obj.bc.(dir)(bcIdx).val{ii} = updateVals{ii};
                end
            end
            
        end
    end
    
    methods % getter functions
        function output = get_boundVal(obj, DIR, tIdx, eIdx, varargin)
            % INPUTS: 
                % varargin -    1st input is rangeMin
                %               2nd input is rangeMax
                %               3rd input is patch depth
            
            % Parse variable inputs for time and element indexing
            % information
            if strcmpi(DIR, 'W') || strcmpi(DIR, 'E')
                nodes = obj.gr.d1vals;
            elseif strcmpi(DIR, 'S') || strcmpi(DIR, 'N')
                nodes = obj.gr.d2vals;
            end
            if ~isempty(varargin)
                if ~isempty(varargin{1})
                    rMin = nodes > varargin{1};
                else
                    rMin = ones(size(nodes));
                end
                if ~isempty(varargin{2})
                    rMax = nodes < varargin{2};
                else
                    rMax = one(size(nodes));
                end
                
                if length(varargin) > 2
                    depth = 2;
                else
                    depth = 1;
                end
            end
            
            outIdx = find(rMin & rMax);
            
            % If time is still empty, select default time
            if isempty(tIdx)
                if length(size(obj.fv)) == 4
                    tIdx = 2;
                else
                    tIdx = 1;
                end
            end
            
            % Perform boundary value extraction
            switch DIR
               case 'W'
                   output = obj.fv(outIdx,depth,:,tIdx);
               case 'E'
                   output = obj.fv(outIdx,end+1 - depth,:,tIdx);
               case 'S'
                   output = obj.fv(depth,outIdx,:,tIdx);
               case 'N'
                   output = obj.fv(end+1 - depth,outIdx,:,tIdx);
            end
            
            % Select element if specified
            if ~isempty(eIdx)
                output = output(:,:,eIdx);
            end
        end
        
        function output = get_boundCondVal(obj, DIR, bcIdx, varargin)
            % OUTPUT: this function outputs the bc values for the given
            % direction
            output = obj.bc.(DIR)(bcIdx).val;
            
            if ~isempty(varargin)
                output = output{varargin{1}};
            end
            
        end
        
    end
    
    methods % helper functions for calculating the gradients
        
        function BC_vals = getBCvals(obj, FF, DIR)
            
            % Setup output to proper size
            if strcmpi(DIR, 'W') || strcmpi(DIR, 'E')
                node_pts = obj.gr.d11(:,1);
            elseif strcmpi(DIR, 'S') || strcmpi(DIR, 'N')
                node_pts = obj.gr.d22(1,:);
            end
            BC_vals = zeros([size(node_pts),size(obj.fv,3)]);
            
            verify = zeros(size(node_pts));
            for i = 1:length(obj.bc.(DIR)) % loop through each subrange
                
                % Extract the corresponding field values
                if ~isempty(obj.bc.(DIR)(i).range)
                    dim_idx = (node_pts > obj.bc.(DIR)(i).range(1)) & (node_pts < obj.bc.(DIR)(i).range(2));
                    if strcmpi(DIR, 'W') || strcmpi(DIR, 'E')
                        FF_sub = FF(dim_idx,:,:);
                    elseif strcmpi(DIR, 'S') || strcmpi(DIR, 'N')
                        FF_sub = FF(:,dim_idx,:);
                    end
                elseif isempty(obj.bc.(DIR)(i).range) || (length(obj.bc.(DIR)) == 1)
                    FF_sub = FF;
                    dim_idx = ones(size(node_pts));
                else
                    error('Ranges not specified!\n');
                end
                
                % calculate partitioned bc
                if strcmp(obj.bc.(DIR)(i).physical, 'sym')
                    BC_temp = obj.gen_symBC(FF_sub, DIR, obj.bc.(DIR)(i).numerical);

                elseif strcmp(obj.bc.(DIR)(i).physical, 'wall')
                    BC_temp = obj.gen_wallBC(FF_sub, DIR, obj.bc.(DIR)(i).numerical, obj.bc.(DIR)(i).val);

                else
                    BC_temp = obj.gen_regBC(FF_sub, DIR, obj.bc.(DIR)(i).numerical, obj.bc.(DIR)(i).val);
                end     
                
                % assign accordingly
                if strcmpi(DIR, 'W') || strcmpi(DIR, 'E')
                    BC_vals(dim_idx,:,:) = BC_temp;
                elseif strcmpi(DIR, 'S') || strcmpi(DIR, 'N')
                    BC_vals(:,dim_idx,:) = BC_temp;
                end
                
                % verify that all points have been assigned?
                verify = verify + dim_idx;
            end
            
            if any(verify > 2)
               error('Redefine ranges! Boundary conditions are overassigned at various points!\nDir: %s \n', DIR);
            elseif any(verify == 0)
                error('Redefine ranges! Boundary conditions are underassigned at various points!\nDir: %s \n', DIR);
            end
            
        end
        
        function BC_vals = gen_symBC(obj, FF, DIR, num_types)
            BC_vals = [];
            for ii = 1:size(FF,3) % loop through vector "elements"
                % symmetry will likely appear right on the
                % boundary of the grid... maybe need to
                % account for displaced symmetry?? 
                switch DIR
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
        
        function BC_vals = gen_regBC(obj, FF, DIR, num_types, vals)
            BC_vals = [];
            for ii = 1:size(obj.fv, 3)
                switch DIR
                    case 'W'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(FF,1),1);
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
                            temp_BC = vals{ii}.*ones(size(FF,1),1);
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
                            temp_BC = vals{ii}.*ones(1, size(FF,2));
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
                            temp_BC = vals{ii}.*ones(1, size(FF,2));
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
        
        function BC_vals = gen_wallBC(obj, FF, DIR, num_types, vals)
            BC_vals = [];
            for ii = 1:size(obj.fv, 3)
                switch DIR
                    case 'W'
                        if isscalar(vals{ii})
                            temp_BC = vals{ii}.*ones(size(FF,1),1);
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
                            temp_BC = vals{ii}.*ones(size(FF,1),1);
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
                            temp_BC = vals{ii}.*ones(1, size(FF,2));
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
                            temp_BC = vals{ii}.*ones(1, size(FF,2));
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