function BC_vals = bcCalc(GR, FL, BC, fv, DIR, varargin)

% BC information stored at each location
% BC.(DIR) contains...
    % - physical type (e.g. wall, inlet, outlet, etc.)
    % - values
    % - variable type (scalar, vector, etc.)
    % - range
    % - dydx (geometry deriative)
    
% BCvals contains the actual boundary values
    % if Dirichlet, then it just returns Dirichlet values
    % if Neumann, evaluates boundary value using linear interp

% Setup output to proper size
%     if strcmpi(DIR, 'W') || strcmpi(DIR, 'E')
%         %node_pts = obj.gr.d11(:,1);
%         BC_vals = zeros(size(fv(:,1,:)));
%     elseif strcmpi(DIR, 'S') || strcmpi(DIR, 'N')
% %         node_pts = obj.gr.d22(1,:);
%         BC_vals = zeros(size(fv(1,:,:)));
%     end

    % Initialize bc_Vals
    switch DIR
        case {'W', 'E'}
            BC_vals = zeros(size(fv,1),1,size(fv,3));
        case {'S', 'N'}
            BC_vals = zeros(1,size(fv,2),size(fv,3));  
    end
    
    % calculate partitioned bc
    for i = 1:length(BC.(DIR))
        if strcmp(BC.(DIR)(i).physical, 'sym')
            BC_temp = gen_symBC(GR, BC, fv, DIR,i);
        elseif strcmp(BC.(DIR)(i).physical, 'wall')
            BC_temp = gen_wallBC(GR, BC, fv, DIR,i);
        elseif strcmp(BC.(DIR)(i).physical, 'outlet')
            BC_temp = gen_outletBC(GR, FL, BC, fv, DIR,i);
        elseif strcmp(BC.(DIR)(i).physical, 'farfield')
            BC_temp = gen_farfieldBC(GR, BC, fv, DIR,i);
        elseif strcmp(BC.(DIR)(i).physical, 'smallDisturb')
            BC_temp = gen_smallDisturbBC(GR, BC, fv, DIR,i);
        else
            BC_temp = gen_inletBC(GR, BC, fv, DIR,i);
        end 
        switch DIR
            case {'W', 'E'}
                BC_vals(BC.(DIR)(i).range(1):BC.(DIR)(i).range(end),:,:) = BC_temp(BC.(DIR)(i).range(1):BC.(DIR)(i).range(end),:,:);
            case {'S', 'N'}
                BC_vals(:,BC.(DIR)(i).range(1):BC.(DIR)(i).range(end),:) = BC_temp(:,BC.(DIR)(i).range(1):BC.(DIR)(i).range(end),:);
        end
    
    end
    
    % first varargin indexes the specific field values interested
    if ~isempty(varargin)
        BC_vals = BC_vals(:,:,varargin{1});
    end


end

function BC_vals = gen_symBC(GR, BC, FF, DIR, ind)
    BC_vals = [];
    % values inputted into BC will represent the axis of symmetry
    switch DIR
        case 'W'
            ind3 = reshape(strcmp(BC.(DIR)(ind).varType, 'v1') + strcmp(BC.(DIR)(ind).varType, 's') + strcmp(BC.(DIR)(ind).varType, 'phi') - strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            if all((GR.XX(:,1)-BC.(DIR)(ind).val)==0) || all((GR.YY(:,1)-BC.(DIR)(ind).val)==0) % -> checks if either YY or XX grid is on axis
                %BC_vals = cat(3, BC_vals, FF(:,2,ii)); 
                BC_vals = FF(:,2,:).*ind3;
            else
                %BC_vals = cat(3, BC_vals, FF(:,1,ii)); % -> XX/YY not on grid axis
                BC_vals = FF(:,1,:).*ind3;
            end

        case 'E'                
            ind3 = reshape(strcmp(BC.(DIR)(ind).varType, 'v1') + strcmp(BC.(DIR)(ind).varType, 's') + strcmp(BC.(DIR)(ind).varType, 'phi') - strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            if all((GR.XX(:,end)-BC.(DIR)(ind).val)==0) || all((GR.YY(:,end)-BC.(DIR)(ind).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(:,end-1,:).*ind3; 
            else
                BC_vals = FF(:,end,:).*ind3; % -> XX/YY not on grid axis
            end 

        case 'N'                
            ind3 = reshape(strcmp(BC.(DIR)(ind).varType, 'v2') + strcmp(BC.(DIR)(ind).varType, 's') - strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            if all((GR.XX(end,:)-BC.(DIR)(ind).val)==0) || all((GR.YY(end,:)-BC.(DIR)(ind).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(end-1,:,:).*ind3; 
            else
                BC_vals = FF(end,:,:).*ind3; % -> XX/YY not on grid axis
            end

        case 'S'
            ind3 = reshape(strcmp(BC.(DIR)(ind).varType, 'v2') + strcmp(BC.(DIR)(ind).varType, 's') - strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            if all((GR.XX(1,:)-BC.(DIR)(ind).val)==0) || all((GR.YY(1,:)-BC.(DIR)(ind).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(2,:,:).*ind3; 
            else
                BC_vals = FF(1,:,:).*ind3; % -> XX/YY not on grid axis
            end 
    end

end

function BC_vals = gen_outletBC(GR, FL, BC, FF, DIR,ind)
    BC_vals = [];
    % values inputted into BC will represent the outlet condition
    % note that the outlet condition is calculated assuming constant spacing
    switch DIR
        case 'W'
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,2,:).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
            BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
%             BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (5./2.*FF(:,1,:) - 2.*FF(:,2,:)+0.5.*FF(:,3,:)).*indPhi;
%             BC_vals = FF(:,1,:);

        case 'E' 
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,end-1,:).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;
%             BC_vals = (4/3.*FF(:,end,:) - 1/3.*FF(:,end-1,:)).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;
            
            if isfield(BC.(DIR), 'exitCond')
                indE = reshape(strcmp(BC.(DIR)(ind).varName, '\rho e'),1,1,size(FF,3));
                indV1 = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
                indV2 = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
                indRho = reshape(strcmp(BC.(DIR)(ind).varName, '\rho'),1,1,size(FF,3));
%                 BC_vals = (5/2.*FF(:,end,:) - 2.*FF(:,end-1,:) + 0.5*FF(:,end-2,:)).*(indTan + indPhi);
                BC_vals = (4/3.*FF(:,end,:) - 1/3.*FF(:,end-1,:)).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;

                if strcmp(BC.(DIR)(ind).exitCond{1}, 'p') % check if exit pressure is defined
                    if any(strcmp(BC.(DIR)(ind).varName, '\rho e')) % full Euler case
                        BC_vals(:,:,indE) = BC.(DIR)(ind).exitCond{2}./(FL.gam-1) + 0.5.*(BC_vals(:,:,indV1).^2 + BC_vals(:,:,indV2).^2)./BC_vals(:,:,indRho);

                    elseif any(strcmp(BC.(DIR)(ind).varName, '\rho u')) % isentropic euler case
%                         indRho = reshape(strcmp(BC.(DIR)(ind).varName, '\rho'),1,1,size(FF,3));
                        BC_vals(:,:,indRho) = ((FL.gam.*FL.M0.*BC.(DIR)(ind).exitCond{2}).^(1/FL.gam));
                    end

                elseif strcmp(BC.(DIR)(ind).exitCond{1}, 'u') % -> define this for exit velocity?
                    if any(strcmp(BC.(DIR)(ind).varName, '\rho e')) % full Euler case
                        BC_vals(:,:,indE) = BC_vals(:,:,indE) - 0.5.*(BC_vals(:,:,indV2).^2)./BC_vals(:,:,indRho) + 0.5.*BC_vals(:,:,indRho).*(BC.(DIR)(ind).exitCond{2}.^2);
                        BC_vals(:,:,indV2) = BC_vals(:,:,indRho).*BC.(DIR)(ind).exitCond{2};
                    end
                end
                
            else
                BC_vals = (4/3.*FF(:,end,:) - 1/3.*FF(:,end-1,:)).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;

            end

        case 'N' 
            BC_vals = FF(end,:,:);
            
        case 'S'
            BC_vals = FF(1,:,:);
            
    end

end

function BC_vals = gen_inletBC(GR, BC, FF, DIR, ind)
    %for ii = 1:size(FF, 3)
    switch DIR
        case {'W','E'}
            BC_vals = zeros(size(FF(:,1,:)));
            
        case {'N', 'S'}
            BC_vals = zeros(size(FF(end,:,:)));
            
    end
    
    for ii = 1:size(FF,3)
        if isscalar(BC.(DIR)(ind).val{ii})
            BC_vals(:,:,ii) = BC.(DIR)(ind).val{ii}; 
        else
            BC_vals(:,:,ii) = reshape(BC.(DIR)(ind).val{ii}, size(BC_vals(:,:,ii))); 
        end
    end

end

function BC_vals = gen_farfieldBC(GR, BC, FF, DIR, ind)
    %for ii = 1:size(FF, 3)
    indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential, expects velocity input instead of actual potential at far-field
    indReal = reshape(~strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3));
    switch DIR
        case 'W'
            BC_vals = FF(:,2,:).*indPhi;
            if GR.isPolar
                dx = GR.RR(:,1).*GR.dT;
            else
                dx = GR.dx;
            end
            
        case 'E'
            BC_vals = FF(:,end-1,:).*indPhi;
            if GR.isPolar
                dx = GR.RR(:,end).*GR.dT;
            else
                dx = GR.dx;
            end
            
        case 'N'
            BC_vals = FF(end-1,:,:).*indPhi;
            if GR.isPolar
                dx = GR.dR;
            else
                dx = GR.dy;
            end
            
        case 'S'
            BC_vals = FF(2,:,:).*indPhi;
            if GR.isPolar
                dx = GR.dR;
            else
                dx = GR.dy;
            end
    end
    
    for ii = 1:size(FF,3)
        if isscalar(BC.(DIR)(ind).val{ii})
            BC_vals(:,:,ii) = BC_vals(:,:,ii) + (indReal(ii) + 2*dx*indPhi(ii))*BC.(DIR)(ind).val{ii}; 
        else
            BC_vals(:,:,ii) = BC_vals(:,:,ii) + (indReal(ii) + 2*dx*indPhi(ii)).*reshape(BC.(DIR)(ind).val{ii}, size(BC_vals(:,:,ii))); 
        end
    end

end

function BC_vals = gen_wallBC(GR, BC, FF, DIR, ind)
    % if curved... assume body is relatively thin
    BC_vals = [];
    switch DIR
        case 'W'
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            BC_vals = 2.*(  FF(:,1,:).*(indTan+indScalar)...
                            - GR.dx.*reshape(BC.(DIR)(ind).dydx,size(FF,1),1).*indPhi...
                            + reshape(BC.(DIR)(ind).dydx,size(FF,1),1).*FF(:,1,indTan).*indNorm) - FF(:,1,:);

        case 'E'                
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            BC_vals = 2.*(  FF(:,end,:).*(indTan+indScalar)...
                            - GR.dx.*reshape(BC.(DIR)(ind).dydx,size(FF,1),1).*indPhi...
                            + reshape(BC.(DIR)(ind).dydx,size(FF,1),1).*FF(:,end,indTan).*indNorm) - FF(:,end,:);

        case 'N'    
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(FF,3));
            BC_vals = (  FF(end,:,:).*(indTan+indScalar+indPhi-indNorm)...
                            - GR.dy.*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indPhi...
                            + 2.*FF(end,:,indRho).*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indNorm);... - FF(1,:,:);

        case 'S'
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(FF,3));
            indNorm = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
            if GR.isPolar
                BC_vals = FF(1,:,:).*(indScalar + indPhi + indTan./(1-0.5.*GR.dR./GR.r_cyl) - indNorm)...
                            + 2.*FF(1,:,indRho).*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indNorm;%.*(~indPhi + indPhi.*GR.RR_S(1,:)./GR.RR_N(1,:))
%                 BC_vals = FF(2,:,:);
            else
                BC_vals = (  FF(1,:,:).*(indTan+indScalar+indPhi-indNorm)...
                            - GR.dy.*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indPhi...
                            + 2.*FF(1,:,indRho).*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indNorm);... - FF(1,:,:);
%                             + 2.*FF(1,:,indTan).*reshape(BC.(DIR)(ind).dydx,1,size(FF,2)).*indNorm);... - FF(1,:,:);
                            
            end
    end
end

function BC_vals = gen_smallDisturbBC(GR,BC, FF, DIR, ind)
    % Small Disturbance Boundary Condition
    % Idea: extrapolate values using outflow, then add perturbation
    BC_walls = [];
    indRho = reshape(strcmp(BC.(DIR)(ind).varName, '\rho'),1,1,size(FF,3));
    
    switch DIR
        case 'W'
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,2,:).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
            BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
%             BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (5./2.*FF(:,1,:) - 2.*FF(:,2,:)+0.5.*FF(:,3,:)).*indPhi;
%             BC_vals = FF(:,1,:);

        case 'E' 
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,end-1,:).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;
            BC_vals = (4/3.*FF(:,end,:) - 1/3.*FF(:,end-1,:)).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;
            
            if isfield(BC.(DIR), 'exitCond')
                if strcmp(BC.(DIR)(ind).exitCond{1}, 'p') % check if exit pressure is defined
                    if any(strcmp(BC.(DIR)(ind).varName, '\rho e')) % full Euler case
                        indE = reshape(strcmp(BC.(DIR)(ind).varName, '\rho e'),1,1,size(FF,3));
                        indV1 = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3));
                        indV2 = reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3));
                        indRho = reshape(strcmp(BC.(DIR)(ind).varName, '\rho'),1,1,size(FF,3));
                        BC_vals(:,:,indE) = BC.(DIR)(ind).exitCond{2}./(FL.gam-1) + 0.5.*(BC_vals(:,:,indV1).^2 + BC_vals(:,:,indV2).^2)./BC_vals(:,:,indRho);

                    elseif any(strcmp(BC.(DIR)(ind).varName, '\rho u')) % isentropic euler case
                        indRho = reshape(strcmp(BC.(DIR)(ind).varName, '\rho'),1,1,size(FF,3));
                        BC_vals(:,:,indRho) = ((FL.gam.*FL.M0.*BC.(DIR)(ind).exitCond{2}).^(1/FL.gam));
                    end

    %             else -> maybe define this for exit density or exit velocity?


                end
            end

        case 'N' 
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,end-1,:).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;
            BC_vals = (4/3.*FF(end,:,:) - 1/3.*FF(end-1,:,:)).*indTan + (2.*FF(end,:,:) - FF(end-1,:,:)).*indPhi;
            
        case 'S'
%             BC_vals = FF(1,:,:);
            indTan = reshape(strcmp(BC.(DIR)(ind).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR)(ind).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR)(ind).varType, 'phi'),1,1,size(FF,3)); % checks for potential
%             BC_vals = FF(:,2,:).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
            BC_vals = (4/3.*FF(1,:,:) - 1/3.*FF(2,:,:)).*indTan + (2.*FF(1,:,:) - FF(2,:,:)).*indPhi;
%             BC_vals = (2.*FF(1,:,:) - FF(2,:,:)).*indTan + (2.*FF(1,:,:) - FF(2,:,:)).*indPhi;
            
    end

    if isfield(BC.(DIR)(ind), 'perturb')
        BC_vals = BC_vals + BC.(DIR)(ind).perturb;
    elseif isfield(BC.(DIR)(ind), 'perturbPrim')
        indices = all(all(BC.(DIR)(ind).perturbPrim ~= 0));
        switch DIR
            case 'E'
                BC_vals = BC_vals + FF(:,end,indRho).*BC.(DIR)(ind).perturbPrim;
            case 'W'
                BC_vals = BC_vals + FF(:,1,indRho).*BC.(DIR)(ind).perturbPrim;
            case 'S'
                BC_vals(:,:,indices) = BC_vals(:,BC.(DIR)(ind).range(1)-1,indices) + FF(1,:,indRho).*BC.(DIR)(ind).perturbPrim(:,:,indices);
            case 'N'
                BC_vals = BC_vals + FF(end,:,indRho).*BC.(DIR)(ind).perturbPrim;
        end
    end

end