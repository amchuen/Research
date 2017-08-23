classdef eulerIsentropicField < handle
% NOTES:
% - This (domain) object manages the field vectors
%   - field vectors can only operate on the elements that make it up, and
%       know info about their own elements,
%   - domain objects manage the vector and have higher level of abstraction
%       with regards to what kind of mathematical variables are stored in the
%       field vectors, as well as their relationships to each other

% - This object solves for the isentropic euler equations
%     - Three-level scheme is the only solving scheme implemented so far
%     - numerical diffusion is introduced to form the shockwave
    
    properties
        GR % grid properties
        FL % fluid properties
        BC % boundary conditions
        CT % control params
        
        % Field Vectors
        FV
        PP % pressure values derived from density
        
    end
    
    methods % methods for creating the object and setting up the field
        function obj = eulerIsentropicField(ctrl, grid_lims, fluid_params, BC_setup)
            % This function helps initialize and build the object
            obj.init_ct(ctrl);
            obj.init_grid(grid_lims); % initialize the grid?
            obj.init_fl(fluid_params); % initialize the 
            obj.init_bc(BC_setup); % initialize boundary conditions?
            obj.init_fv;
        end 
        
        % functions that help smooth operation of the class, not meant to
        % be used by user
        
        function obj = init_ct(obj, ctrl)
            
            % numerical controls
            obj.CT.eps_s = ctrl.eps_s; % spatial diffusion term
            obj.CT.eps_t = ctrl.eps_t; % time diffusion term
            obj.CT.tol = ctrl.tol;
            obj.CT.dt = ctrl.dt;
            obj.CT.iter_min = ctrl.iter_min;
            obj.CT.CFL_on = ctrl.CFL_on;
            obj.CT.use_1visc = ctrl.use_1visc;
            obj.CT.is_polar = ctrl.is_polar;
            obj.CT.case_name = ctrl.case_name;
            
            % object controls
            
        end
        
        function obj = init_grid(obj, grid_lims)
            obj.GR.d1 = grid_lims(1,3);
            obj.GR.d2 = grid_lims(2,3);
            obj.GR.d1vals = grid_lims(1,1):obj.GR.d1:grid_lims(1,2);
            obj.GR.d2vals = grid_lims(2,1):obj.GR.d2:grid_lims(2,2);
            [obj.GR.d22, obj.GR.d11] = meshgrid(obj.GR.d2vals, obj.GR.d1vals);
            
            if obj.CT.is_polar
                obj.GR.XX = obj.GR.d22 .* cos(obj.GR.d11);
                obj.GR.YY = obj.GR.d22 .* sin(obj.GR.d11);
            else
                obj.GR.XX = obj.GR.d22;
                obj.GR.YY = obj.GR.d11;
            end
            
        end
        
        function obj = init_fl(obj, fluid_params)
           obj.FL.M0 = fluid_params.M0;
           obj.FL.gam = fluid_params.gam;
        end
        
        function obj = init_bc(obj, BC_setup)
            % assumes BC_setup is a cell value, with first row the grid
            % face name, the second row the BC type, and the third row the
            % BC values in vector form for rho, rho*u, and rho*v
            % Different types of BC's
            % * Wall -> N or D
            % * Inlet -> D
            % * Outflow -> N
            % * Far-field -> N
            % * Symmetry/Reflection -> N or D
            % * Patch between fields -> D            
            
%             obj.BC = BC_setup;
            for vec = 1:3
                for i = 1:size(BC_setup,2)
                    % assign the BC type to the struct title
                    % sets up BC type for each fieldVector object
                    % what is needed: direction, value, BC type
                    % value implied from BC type?
                    
                    obj.BC(vec).(BC_setup{1,i}).physical = BC_setup{2,i};
                    BC_vals = {};
                    switch BC_setup{2,i} % will need to modify to account for polar cases
                        case 'wall'
                            
                            % determine which velocity is normal to the
                            % wall
                            if (strcmp(BC_setup{1,i},'W') || strcmp(BC_setup{1,i}, 'E'))
                                wall_dir = 2;
                            elseif (strcmp(BC_setup{1,i},'S') || strcmp(BC_setup{1,i}, 'N'))
                                wall_dir = 3;
                            end
                            
                            % check if wall is zero or a linearized BC
                            if all(BC_setup{wall_dir+2,i} == 0)
                                obj.BC(vec).(BC_setup{1,i}).update = 0;
                            else
                                obj.BC(vec).(BC_setup{1,i}).update = 1; % update must be done if wall has a non-zero shape
                            end
                            
                            if vec == wall_dir
                                obj.BC(vec).(BC_setup{1,i}).numerical = {'D', 'D', 'D'};
                                obj.BC(vec).(BC_setup{1,i}).val = {BC_setup{vec+2,i}, BC_setup{vec+2,i}, BC_setup{vec+2,i}};
                                obj.BC(vec).(BC_setup{1,i}).val{vec} = obj.BC(vec).(BC_setup{1,i}).val{vec}.*BC_setup{vec+2,i};
                                
                            else
                                obj.BC(vec).(BC_setup{1,i}).numerical = {'N', 'N', 'N'};
                                obj.BC(vec).(BC_setup{1,i}).numerical{wall_dir} = 'D';
                                obj.BC(vec).(BC_setup{1,i}).val = BC_setup(3:5,i)';
                            end
                            
                        case 'inlet'
                            obj.BC(vec).(BC_setup{1,i}).update = 0;
                            obj.BC(vec).(BC_setup{1,i}).numerical = {'D', 'D', 'D'};
                            
                            for ii = 3:size(BC_setup,1)
                               BC_vals = cat(2, BC_vals, BC_setup{ii,i}.*BC_setup{vec+2,i}./BC_setup{3,i});
                            end
                            obj.BC(vec).(BC_setup{1,i}).val = BC_vals;
                            
                        case 'outlet'
                            obj.BC(vec).(BC_setup{1,i}).update = 0;
                            obj.BC(vec).(BC_setup{1,i}).numerical = {'N', 'N', 'N'};
                            
                            obj.BC(vec).(BC_setup{1,i}).val = {0,0,0};
                            
                        case 'far-field'
                            obj.BC(vec).(BC_setup{1,i}).update = 0;
                            obj.BC(vec).(BC_setup{1,i}).numerical = {'D', 'D', 'D'};
                            
                            for ii = 3:size(BC_setup,1)
                                BC_vals = cat(2, BC_vals, BC_setup{ii,i}.*BC_setup{vec+2,i}./BC_setup{3,i});
                            end
                            obj.BC(vec).(BC_setup{1,i}).val = BC_vals;
                            
                        case 'sym'
                            % in this case, if the numerical type is
                            % Neumann, then it is implied that the value is
                            % not normal to the surface, wherewas Dirichlet
                            % will be; i.e. the two distinctions will serve
                            % a different purpose
                            
                            % need to check if boundary is on axis or not?
                            obj.BC(vec).(BC_setup{1,i}).update = 0;
                            % updates not necessary, as can be extrapolated
                            % from the field values
                            if (strcmp(BC_setup{1,i},'W') || strcmp(BC_setup{1,i}, 'E'))
                                obj.BC(vec).(BC_setup{1,i}).numerical = {'N', 'D', 'N'};
                            elseif (strcmp(BC_setup{1,i},'S') || strcmp(BC_setup{1,i}, 'N'))
                                obj.BC(vec).(BC_setup{1,i}).numerical = {'N', 'N', 'D'};
                            end
                            
                        case 'patch'
                            obj.BC(vec).(BC_setup{1,i}).update = 1;
                            obj.BC(vec).(BC_setup{1,i}).numerical = {'D', 'D', 'D'};
                            
                            for ii = 3:size(BC_setup,1)
                               BC_vals = cat(2, BC_vals, BC_setup{ii,i}.*BC_setup{vec+2,i}./BC_setup{3,i});
                            end
                            obj.BC(vec).(BC_setup{1,i}).val = BC_vals;
                            
                    end
                    
                end
                
            end
            
        end
        
        function obj = init_fv(obj)
           % initialize vectorized variables for solving isentropic euler 
            fv = {ones(size(obj.GR.d22)), ones(size(obj.GR.d22)), zeros(size(obj.GR.d22))};
            for i = 1:3
                if i == 1
                    obj.FV = fieldVector(obj.GR, obj.BC(i), obj.CT, 1, fv); 
                else
                    fv_i = {fv{i}, fv{2}.*fv{i}./fv{1}, fv{3}.*fv{i}./fv{1}};
                    obj.FV(i) = fieldVector(obj.GR, obj.BC(i), obj.CT, 0, fv_i); 
                end
            end
            
            % initialize residual array
            obj.CT.res = [];%zeros(size(obj.FV));

            PP_fv = obj.FV(1).fv(:,:,1,2).^(obj.FL.gam) ./ (obj.FL.gam .* obj.FL.M0.^2); % pressure 
            % need to update bc conditions for pressure here!
            fnames = fieldnames(obj.BC(1));
            obj.BC(4) = obj.BC(1);
            for ii = 1:length(fnames)
                obj.BC(4).(fnames{ii}).val{1} = obj.BC(4).(fnames{ii}).val{1} .^ (obj.FL.gam) ./ (obj.FL.gam .* obj.FL.M0.^2);
            end
            
            obj.PP = fieldScalar(obj.GR, obj.BC(4), obj.CT, PP_fv);
        end
    end
    
    methods % solving schemes
        
        function obj = timeStep_tl(obj)

            % Calculate CFL condition and other stability parameters
            obj = obj.CFLcheck;

            % Perform time-step
            for i = 1:length(obj)
                P_grad = cat(3, zeros(size(obj(i).PP.fv)), obj(i).PP.f_2, obj(i).PP.f_1); % pressure gradients for momentum eqns

        %                 if obj(i).CT.is_polar
        %                     alpha = cat(3, zeros(size(obj(i).GR.RR)), 1./(obj(i).GR.RR(:,:,2:3).^2));
        %                     DF_coeff = ((rr_n + rr_s)./(dr^2 .* RR) + 2./(RR.^2 .* dT^2) + alpha);
        %                 else
        %                 DF_coeff =  (2./(obj(i).GR.d1^2) + 2./(obj(i).GR.d2^2));
        %                 end

                if obj(i).CT.use_1visc
                    DF_coeff =  (2./(obj(i).GR.d1^2) + 2./(obj(i).GR.d2^2));
                    obj(i).FV(1).fv(:,:,:,3) = (obj(i).CT.eps_s.*(obj(i).FV(1).laplace + DF_coeff.*obj(i).FV(1).fv(:,:,:,2)) - 0.5.*obj(i).FV(1).fv(:,:,:,1).*(obj(i).CT.eps_s.*DF_coeff - 1./obj(i).CT.dt) - P_grad - (obj(i).FV(2).f_2 + obj(i).FV(3).f_1))./...
                                (0.5./obj(i).CT.dt + 0.5.*obj(i).CT.eps_s.*DF_coeff);
                else
                    obj(i).FV(1).fv(:,:,:,3) = (obj(i).CT.eps_s .* obj(i).FV(1).laplace - P_grad - (obj(i).FV(2).f_2 + obj(i).FV(3).f_1) +...
                                            2.*obj(i).CT.eps_t.*obj(i).FV(1).fv(:,:,:,2)./(obj(i).CT.dt^2) - (obj(i).CT.eps_t./(obj(i).CT.dt.^2) - 0.5./(obj(i).CT.dt)).*obj(i).FV(1).fv(:,:,:,1))./...
                                            (obj(i).CT.eps_t./(obj(i).CT.dt^2) + 0.5./(obj(i).CT.dt));
                end

                obj(i).fvChecks();

            end
            
            % Output residual if desired
            obj = obj.calcRes(0);
            
        end
        
        function fvChecks(obj)    
            % Perform checks
            if any(any(max(abs(obj.FV(1).fv(:,:,3,end)./obj.FV(1).fv(:,:,1,end)))>10))
               fprintf('Check Tangential veocity!\n');
            end

            if (~isreal(obj.FV(1).fv(:,:,:,end)) || any(any(any(isnan(obj.FV(1).fv(:,:,:,end))))))
                fprintf('Solution exhibits non-solutions (either non-real or NaN) in nodes!\n');
%                 break;
            end
        end
        
        function out = checkConvergence(obj)
            out = zeros(size(obj));
            for i = 1:length(obj)
                out(i) = isempty(obj(i).CT.res) || ((max(obj(i).CT.res(end, :)) > obj(i).CT.tol*max(obj(i).CT.res(obj(i).CT.res<1)))|| (size(obj(i).CT.res,1) < obj(i).CT.iter_min));
            end
            out = any(out);
        end
        
        function obj = updateVals(obj)
            for ii = 1:length(obj)
                obj(ii) = obj(ii).updateFV_tl;
            end
            
            for ii = 1:length(obj)
                obj(ii) = obj(ii).updateBC_tl(obj((1:length(obj)) ~= ii));

                for i = 1:length(obj(ii).FV)
                    obj(ii).FV(i) = obj(ii).FV(i).Gradient();
                    obj(ii).FV(i) = obj(ii).FV(i).Laplacian();
                end
            end               
        end
        
        function obj = updateFV_tl(obj)
            % updates the object
            % update procedure should include methods to update the derivs
            % once values are updated
            for i = 1:length(obj.FV) % loop through vectors
                if i == 1
                   obj.FV(i) = obj.FV(i).update_fv();
                else
                   argin = cell(1,size(obj.FV(i).fv,3));
                   for ii = 1:size(obj.FV(i).fv,3) % loop through elems
                      argin{ii} = obj.FV(1).fv(:,:,ii,2) .* obj.FV(1).fv(:,:,i,2)./obj.FV(1).fv(:,:,1,2); 
                   end
                   obj.FV(i) = obj.FV(i).update_fv(argin);
                   
                end                
            end
            
            obj.PP = obj.PP.update_fv(obj.FV(1).fv(:,:,1,2) .^ obj.FL.gam ./ (obj.FL.gam .* obj.FL.M0^2));
            
        end
        
        function obj = updateBC_tl(obj, others)
            % NOTES:
            % - this function is run after field values have been updated
            % - densities should always initialized with 1!
            
            fnames = fieldnames(obj.BC);
            for i = 1:length(fnames) % loop through each direction
                
                if obj.BC(1).(fnames{i}).update
                    
                    switch obj.BC(1).(fnames{i}).physical
                        case 'wall'
                            obj = obj.updateWallBC(fnames{i});
                            
                        case 'patch'
                            
                            % first match the boundaries
                            % assume each cell of varargin represents the
                            % other corresponding objects
                            
                            % test for matrix size, array distance, and
                            % opposite side
                            for fv = 1:length(others)
                                switch fnames{i}
                                    case 'W'
                                        test_dir = 'E';
                                        test_home = obj.GR.d22(:,1);
                                        test_patch = others(fv).GR.d22(:,end);
                                        test_dx = obj.GR.d2;
                                    case 'E'
                                        test_dir = 'W';
                                        test_home = obj.GR.d22(:,end);
                                        test_patch = others(fv).GR.d22(:,1);
                                        test_dx = obj.GR.d2;
                                    case 'S'
                                        test_dir = 'N';
                                        test_home = obj.GR.d11(1,:);
                                        test_patch = others(fv).GR.d11(end,:);
                                        test_dx = obj.GR.d1;
                                    case 'N'
                                        test_dir = 'S';
                                        test_home = obj.GR.d11(end,:);
                                        test_patch = others(fv).GR.d11(1,:);
                                        test_dx = obj.GR.d1;
                                end
                                
                                if (length(test_home) ~= length(test_patch))
                                    continue;
                                    
                                elseif (abs(max(abs(test_home - test_patch)) - test_dx) < 1e-8)
%                                     BC_vals = others(fv).FV(1).get_boundVal(test_dir);
                                    break;
                                end
                                
                            end
                            
                            % update boundary condition
                            for vec = 1:length(obj.FV) % loop through vectors
                                BC_vals = others(fv).FV(vec).get_boundVal(test_dir);
                                obj.FV(vec) = obj.FV(vec).update_bc(fnames{i}, mat2cell(BC_vals, size(BC_vals,1), size(BC_vals,2), ones(1, size(BC_vals,3))));
                            end                           
                    end
                end
            end            
        end
        
        function obj = updateWallBC(obj, dir)
            % must update rho * vel for normal velocity
            % then apply to other vectors
            
            % get old density
            d0 = obj.FV(1).get_boundVal(dir, 1,1); % t1

            % get wall normal velocity for updating
            wall_loc = strcmpi(obj.BC(1).(dir).numerical, 'D');
            wall_dir = find(wall_loc);
            for vec = 1:length(obj.FV)
                update_vals = [];
                update_cell = cell(1,length(obj.FV));    
                
                if vec ~= wall_dir
                    % update for element direction that is normal to the wall
                    update_vals = obj.FV(1).get_boundVal(dir,[],vec).*obj.FV(1).get_boundCondVal(dir, wall_dir)./d0; 
                    update_cell{wall_dir} = update_vals;
                else
                    % we are at the vector that corresponds to normal wall vel
                    update_vals = obj.FV(1).get_boundVal(dir).*repmat(obj.FV(1).get_boundCondVal(dir, wall_dir)./d0, 1,1,size(obj.FV(1).fv,3)); 
                    update_cell = mat2cell(update_vals, size(update_vals,1), size(update_vals,2), ones(1,size(update_vals,3)));
                end
                
                obj.FV(vec) = obj.FV(vec).update_bc(dir, update_cell);
            end
        end
        
        function obj = CFLcheck(obj)
            CFL_i = 0;
            for i = 1:length(obj)
               % Check CFL conditions
                Ur = (obj(i).FV(1).fv(:,:,2,2)./obj(i).FV(1).fv(:,:,1,2));
                VT = (obj(i).FV(1).fv(:,:,3,2)./obj(i).FV(1).fv(:,:,1,2));
                CFL = (max(abs(Ur(:)))./(obj(i).GR.d2) + max(abs(VT(:)))./obj(i).GR.d1)*obj(i).CT.dt;
                if CFL > CFL_i % get maximum CFL value for entire simulation
                    CFL_i = CFL;
                end
            end

            if CFL_i >= 1.0 %&& CFL_on
               fprintf('CFL condition not met!\n');
               fprintf('Decreasing time steps!\n');
               for i = 1:length(obj)
                    obj(i).CT.dt = obj(i).CT.dt*0.8 / CFL_i;
               end
               fprintf('New time step:%0.5f\n', obj(1).CT.dt);
            end 
        end
          
        function obj = calcRes(obj, plotTrigger)
            for i = 1:length(obj)
                err = abs(obj(i).FV(1).fv(:,:,:,3) - obj(i).FV(1).fv(:,:,:,2));

                if isempty(obj(i).CT.res) %(size(obj(i).CT.res,1) == 1) && all(obj(i).CT.res(end,:) == 0)
                    for ii = 1:size(obj(i).FV,2)
                        obj(i).CT.res(1, ii) = max(max(err(:,:,ii))); 
                    end
                else
                    obj(i).CT.res(end+1, :) = zeros(size(obj(i).FV));
                    for ii = 1:size(obj(i).FV,2)
                        obj(i).CT.res(end, ii) = max(max(err(:,:,ii))); 
                    end

    %                 obj.CT.res(end+1, :) = [max(R_err(:)), max(A_err(:)), max(B_err(:))];
                end
            end
            

            if ((size(obj(1).CT.res, 1) > 500) && (mod(size(obj(1).CT.res, 1), 2000) == 0)) || (plotTrigger)
                whole_res = [];
                for i = 1:length(obj)
                   whole_res = cat(3, whole_res, obj(i).CT.res); 
                end
                
                whole_res = max(whole_res,[], 3);
                
                fprintf('Iteration Ct: %i\n', size(whole_res, 1));
                fprintf('Current Residual: %0.5e\n', max(whole_res(end, :)));
                toc;
                figure(1);semilogy(1:size(whole_res,1), whole_res(:,1));
                hold on;
                semilogy(1:size(whole_res,1), whole_res(:,2));
                semilogy(1:size(whole_res,1), whole_res(:,3));
                hold off;
                legend('Density', '\rho u', '\rho v', 'Location' ,'bestoutside');
                fprintf('\n');
            end
        end
        
    end
    
    methods % methods for plotting and displaying results
        function post_process(obj, DIR)

%             Ux = (obj.FV(1).fv(:,:,2,3) ./ obj.FV(1).fv(:,:,1,3));
%             Vy = (obj.FV(1).fv(:,:,3,3) ./ obj.FV(1).fv(:,:,1,3));
% 
%             q2_ij = (Ux).^2 + (Vy).^2;
            
            whole_res = [];
            for i = 1:length(obj)
                whole_res = cat(3, whole_res, obj(i).CT.res);
            end
            
            whole_res = max(whole_res,[], 3);

            close all;
            figure(1);
            for i = 1:size(whole_res, 2)
            semilogy(1:size(whole_res,1), whole_res(:,i));
            hold on;
            end
            hold off;
            
            legend('Density', '\rho u', '\rho v');
            title(['Residual Plot, M=' num2str(obj(1).FL.M0)]);
            xlabel('# of iterations');
            ylabel('Residual (Error)');

            saveas(gcf, [DIR '\residual_plot.pdf']);
            saveas(gcf, [DIR '\residual_plot']);
%             
%              figure(2);
%             % plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
%             plot(obj(1).GR.XX(1,:), 1-q2_ij(1,:));
            
            %% contour plots
%             XX1 = cat(1, (obj(1).GR.XX), (obj(2).GR.XX));
%             XX2 = cat(2, (obj(2).GR.XX(:,end)), (obj(3).GR.XX));
%             YY1 = cat(1, (obj(1).GR.YY), (obj(2).GR.YY));
%             YY2 = cat(2, (obj(2).GR.YY(:,end)), (obj(3).GR.YY));
%             
%             RHO1 = cat(1, (obj(1).FV(1).fv(:,:,1,3)), (obj(2).FV(1).fv(:,:,1,3)));
%             RHO2 = cat(2, (obj(2).FV(1).fv(:,end,1,3)), (obj(3).FV(1).fv(:,:,1,3)));
%             Ux1 = cat(1, (obj(1).FV(1).fv(:,:,2,3) ./ obj(1).FV(1).fv(:,:,1,3)), (obj(2).FV(1).fv(:,:,2,3) ./ obj(2).FV(1).fv(:,:,1,3)));
%             Vy1 = cat(1, (obj(1).FV(1).fv(:,:,3,3) ./ obj(1).FV(1).fv(:,:,1,3)), (obj(2).FV(1).fv(:,:,3,3) ./ obj(2).FV(1).fv(:,:,1,3)));
%             Ux2 = cat(2, (obj(2).FV(1).fv(:,end,2,3) ./ obj(2).FV(1).fv(:,end,1,3)), (obj(3).FV(1).fv(:,:,2,3) ./ obj(3).FV(1).fv(:,:,1,3)));
%             Vy2 = cat(2, (obj(2).FV(1).fv(:,end,3,3) ./ obj(2).FV(1).fv(:,end,1,3)), (obj(3).FV(1).fv(:,:,3,3) ./ obj(3).FV(1).fv(:,:,1,3)));
%             P1 = cat(1, (obj(1).PP.fv), (obj(2).PP.fv));
%             P2 = cat(2, (obj(2).PP.fv(:,end)), (obj(3).PP.fv));
%             
%             q2_ij1 = (Ux1).^2 + (Vy1).^2;
%             q2_ij2 = (Ux2).^2 + (Vy2).^2;

%             plot density
            for i = 1:length(obj)
                Ux = (obj(i).FV(1).fv(:,:,2,3) ./ obj(i).FV(1).fv(:,:,1,3));
                Vy = (obj(i).FV(1).fv(:,:,3,3) ./ obj(i).FV(1).fv(:,:,1,3));

                q2_ij = (Ux).^2 + (Vy).^2;
                
                figure(2);
                % plot([fliplr(XX(:,end)'), XX(1,:)], 1 - [fliplr(q2_ij(:,end)'), q2_ij(1,:)]);
                plot(obj(1).GR.XX(1,:), 1-q2_ij(1,:));
                hold on;
                
                figure(3);
                hold on;
                contourf(obj(i).GR.XX,obj(i).GR.YY,round(obj(i).FV(1).fv(:,:,1,2),3), 50)
                

                figure(4); % cp plots
                hold on;
                contourf(obj(i).GR.XX, obj(i).GR.YY, 1-q2_ij, 50); %./((RR.*cos(TT)).^2)

                figure(5); % pressure
                hold on;
                contourf(obj(i).GR.XX, obj(i).GR.YY, obj(i).PP.fv, 50);

                figure(6); % theta-dir velocity plots
                hold on;
                contourf(obj(i).GR.XX, obj(i).GR.YY, Ux, 50); %./((RR.*cos(TT)).^2)
                

                figure(7); % theta-dir velocity plots
                hold on;
                contourf(obj(i).GR.XX, obj(i).GR.YY, Vy, 50); %./((RR.*cos(TT)).^2)
            
            end
            
            figure(2);
            xlabel('X');
            %     ylabel('\phi_{\theta}');
            ylabel('C_p');
            title(['C_p on surface, M=' num2str(obj.FL.M0)]);
            set(gca, 'Ydir', 'reverse');
            saveas(gcf, [DIR '\cp_surf.pdf']);
            saveas(gcf, [DIR '\cp_surf']);
            
            figure(3);
%             contourf(XX1, YY1, RHO1, 50);hold on;
%             contourf(XX2, YY2, RHO2, 50);
            title(['Density (Normalized), M=' num2str(obj(i).FL.M0)]);
            colorbar('eastoutside');
            axis equal
            saveas(gcf, [DIR '\density.pdf']);
            saveas(gcf, [DIR '\density']);
            
            figure(4);
%             contourf(XX1, YY1, q2_ij1, 50);hold on;
%             contourf(XX2, YY2, q2_ij2, 50);
            title(['Pressure Coefficient Contours, M=' num2str(obj(i).FL.M0)]);
            colorbar('eastoutside');
            axis equal
            saveas(gcf, [DIR '\cp_contour.pdf']);
            saveas(gcf, [DIR '\cp_contour']);
            
            figure(5);
%             contourf(XX1, YY1, P1, 50);hold on;
%             contourf(XX2, YY2, P2, 50);
            title(['Pressure (Normalized), M=' num2str(obj(i).FL.M0)]);
            colorbar('eastoutside');
            axis equal
            saveas(gcf, [DIR '\pressure.pdf']);
            saveas(gcf, [DIR '\pressure']);
            
            figure(6);
%             contourf(XX1, YY1, Ux1, 50);hold on;
%             contourf(XX2, YY2, Ux2, 50);
            title(['U velocity, M=' num2str(obj(i).FL.M0)]);
            colorbar('eastoutside');
            axis equal
            saveas(gcf, [DIR '\phi_theta.pdf']);
            saveas(gcf, [DIR '\phi_theta']);
            
            figure(7);
%             contourf(XX1, YY1, Vy1, 50);hold on;
%             contourf(XX2, YY2, Vy2, 50);
            title(['V velocity, M=' num2str(obj(i).FL.M0)]);
            colorbar('eastoutside');
            axis equal
            saveas(gcf, [DIR '\phi_radius.pdf']);
            saveas(gcf, [DIR '\phi_radius']);

            % Save Results
            save([DIR '\results.mat']);
            
        end
        
        function plot_fv(obj, vecInd, elem)
            
            titles = {'\rho', 'A', 'B'};
            
            if nargin == 2
                if vecInd ~= 1
                    for i = 1:size(obj.FV(vecInd).fv,3)
                       figure();
                       contourf(obj.GR.XX, obj.GR.YY, obj.FV(vecInd).fv(:,:,i), 50);
                       if i == 1
                            title(titles{vecInd});
                       else
                            title([titles{vecInd} '*' titles{i} '/' titles{1}]);
                       end
                       colorbar('eastoutside');
                    end
                   
                else
                    for i = 1:size(obj.FV(vecInd).fv,3)
                        figure();
                        contourf(obj.GR.XX, obj.GR.YY, obj.FV(vecInd).fv(:,:,i,3), 50);
                        title(titles{i});
                        colorbar('eastoutside');
                    end
                    
                end
            elseif nargin == 3
                figure();
                if vecInd == 1
                    contourf(obj.GR.XX, obj.GR.YY, obj.FV(vecInd).fv(:,:,elem,3), 50); 
                else
                    contourf(obj.GR.XX, obj.GR.YY, obj.FV(vecInd).fv(:,:,elem), 50);
                end
                if elem == 1
                    title(titles{vecInd});
                else
                    title([titles{vecInd} '*' titles{elem} '/' titles{1}]);
                end
                colorbar('eastoutside');
            end
        end
%         
%         function plot_derivs(vecInd, elem, dir)
%             
%         end
        
    end
end