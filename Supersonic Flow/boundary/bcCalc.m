function BC_vals = bcCalc(GR, BC, fv, DIR, varargin)

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

    % calculate partitioned bc
    
    if strcmp(BC.(DIR).physical, 'sym')
        BC_temp = gen_symBC(GR, BC, fv, DIR);
    elseif strcmp(BC.(DIR).physical, 'wall')
        BC_temp = gen_wallBC(GR, BC, fv, DIR);
    elseif strcmp(BC.(DIR).physical, 'outlet')
        BC_temp = gen_outletBC(GR, BC, fv, DIR);
    elseif strcmp(BC.(DIR).physical, 'farfield')
        BC_temp = gen_farfieldBC(GR, BC, fv, DIR);
    else
        BC_temp = gen_inletBC(GR, BC, fv, DIR);
    end 
    
    BC_vals = BC_temp;
    
    % first varargin indexes the specific field values interested
    if ~isempty(varargin)
        BC_vals = BC_vals(:,:,varargin{1});
    end


end

function BC_vals = gen_symBC(GR, BC, FF, DIR)
    BC_vals = [];
    % values inputted into BC will represent the axis of symmetry
    switch DIR
        case 'W'
            ind3 = reshape(strcmp(BC.(DIR).varType, 'v1') + strcmp(BC.(DIR).varType, 's') + strcmp(BC.(DIR).varType, 'phi') - strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            if all((GR.XX(:,1)-BC.(DIR).val)==0) || all((GR.YY(:,1)-BC.(DIR).val)==0) % -> checks if either YY or XX grid is on axis
                %BC_vals = cat(3, BC_vals, FF(:,2,ii)); 
                BC_vals = FF(:,2,:).*ind3;
            else
                %BC_vals = cat(3, BC_vals, FF(:,1,ii)); % -> XX/YY not on grid axis
                BC_vals = FF(:,1,:).*ind3;
            end

        case 'E'                
            ind3 = reshape(strcmp(BC.(DIR).varType, 'v1') + strcmp(BC.(DIR).varType, 's') + strcmp(BC.(DIR).varType, 'phi') - strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            if all((GR.XX(:,end)-BC.(DIR).val)==0) || all((GR.YY(:,end)-BC.(DIR).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(:,end-1,:).*ind3; 
            else
                BC_vals = FF(:,end,:).*ind3; % -> XX/YY not on grid axis
            end 

        case 'N'                
            ind3 = reshape(strcmp(BC.(DIR).varType, 'v2') + strcmp(BC.(DIR).varType, 's') - strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            if all((GR.XX(end,:)-BC.(DIR).val)==0) || all((GR.YY(end,:)-BC.(DIR).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(end-1,:,:).*ind3; 
            else
                BC_vals = FF(end,:,:).*ind3; % -> XX/YY not on grid axis
            end

        case 'S'
            ind3 = reshape(strcmp(BC.(DIR).varType, 'v2') + strcmp(BC.(DIR).varType, 's') - strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            if all((GR.XX(1,:)-BC.(DIR).val)==0) || all((GR.YY(1,:)-BC.(DIR).val)==0) % -> checks if either YY or XX grid is on axis
                BC_vals = FF(2,:,:).*ind3; 
            else
                BC_vals = FF(1,:,:).*ind3; % -> XX/YY not on grid axis
            end 
    end

end

function BC_vals = gen_outletBC(GR, BC, FF, DIR)
    BC_vals = [];
    % values inputted into BC will represent the outlet condition
    % note that the outlet condition is calculated assuming constant spacing
%     switch DIR
%         case 'W'
%             BC_vals = FF(:,2,:);
% 
%         case 'E'                
%             BC_vals = FF(:,end-1,:);
% 
%         case 'N' 
%             BC_vals = FF(end-1,:,:);
%             
%         case 'S'
%             BC_vals = FF(2,:,:);
%             
%     end

    switch DIR
        case 'W'
            indTan = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            BC_vals = FF(:,1,:).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
%             BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (2.*FF(:,1,:) - FF(:,2,:)).*indPhi;
%             BC_vals = (4/3.*FF(:,1,:) - 1/3.*FF(:,2,:)).*indTan + (5./2.*FF(:,1,:) - 2.*FF(:,2,:)+0.5.*FF(:,3,:)).*indPhi;
%             BC_vals = FF(:,1,:);

        case 'E' 
            indTan = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3)) + reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            BC_vals = FF(:,end-1,:).*indTan + (2.*FF(:,end,:) - FF(:,end-1,:)).*indPhi;

        case 'N' 
            BC_vals = FF(end,:,:);
            
        case 'S'
            BC_vals = FF(1,:,:);
            
    end

end

function BC_vals = gen_inletBC(GR, BC, FF, DIR)
    %for ii = 1:size(FF, 3)
    switch DIR
        case {'W','E'}
            BC_vals = zeros(size(FF(:,1,:)));
            
        case {'N', 'S'}
            BC_vals = zeros(size(FF(end,:,:)));
            
    end
    
    for ii = 1:size(FF,3)
        if isscalar(BC.(DIR).val{ii})
            BC_vals(:,:,ii) = BC.(DIR).val{ii}; 
        else
            BC_vals(:,:,ii) = reshape(BC.(DIR).val{ii}, size(BC_vals(:,:,ii))); 
        end
    end

end

function BC_vals = gen_farfieldBC(GR, BC, FF, DIR)
    %for ii = 1:size(FF, 3)
    indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential, expects velocity input instead of actual potential at far-field
    indReal = reshape(~strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3));
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
        if isscalar(BC.(DIR).val{ii})
            BC_vals(:,:,ii) = BC_vals(:,:,ii) + (indReal(ii) + 2*dx*indPhi(ii))*BC.(DIR).val{ii}; 
        else
            BC_vals(:,:,ii) = BC_vals(:,:,ii) + (indReal(ii) + 2*dx*indPhi(ii)).*reshape(BC.(DIR).val{ii}, size(BC_vals(:,:,ii))); 
        end
    end

end

function BC_vals = gen_wallBC(GR, BC, FF, DIR)
    % if curved... assume body is relatively thin
    BC_vals = [];
    switch DIR
        case 'W'
            indTan = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            BC_vals = 2.*(  FF(:,1,:).*(indTan+indScalar)...
                            - GR.dx.*reshape(BC.(DIR).dydx,size(FF,1),1).*indPhi...
                            + reshape(BC.(DIR).dydx,size(FF,1),1).*FF(:,1,indTan).*indNorm) - FF(:,1,:);

        case 'E'                
            indTan = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            BC_vals = 2.*(  FF(:,end,:).*(indTan+indScalar)...
                            - GR.dx.*reshape(BC.(DIR).dydx,size(FF,1),1).*indPhi...
                            + reshape(BC.(DIR).dydx,size(FF,1),1).*FF(:,end,indTan).*indNorm) - FF(:,end,:);

        case 'N'    
            indTan = reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indNorm = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            BC_vals = 2.*(  FF(end,:,:).*(indTan+indScalar)...
                            - GR.dy.*reshape(BC.(DIR).dydx,1,size(FF,2)).*indPhi...
                            + reshape(BC.(DIR).dydx,1,size(FF,2)).*FF(end,:,indTan).*indNorm) - FF(end,:,:);

        case 'S'
            indTan = reshape(strcmp(BC.(DIR).varType, 'v2'),1,1,size(FF,3));
            indScalar = reshape(strcmp(BC.(DIR).varType, 's'),1,1,size(FF,3));
            indPhi = reshape(strcmp(BC.(DIR).varType, 'phi'),1,1,size(FF,3)); % checks for potential
            indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(FF,3));
            indNorm = reshape(strcmp(BC.(DIR).varType, 'v1'),1,1,size(FF,3));
            if GR.isPolar
                BC_vals = FF(1,:,:).*(indScalar + indPhi + indTan./(1-0.5.*GR.dR./GR.r_cyl) - indNorm)...
                            + 2.*FF(1,:,indScalar).*reshape(BC.(DIR).dydx,1,size(FF,2)).*indNorm;%.*(~indPhi + indPhi.*GR.RR_S(1,:)./GR.RR_N(1,:))
%                 BC_vals = FF(2,:,:);
            else
                BC_vals = (  FF(1,:,:).*(indTan+indScalar+indPhi-indNorm)...
                            - GR.dy.*reshape(BC.(DIR).dydx,1,size(FF,2)).*indPhi...
                            + 2.*FF(1,:,indRho).*reshape(BC.(DIR).dydx,1,size(FF,2)).*indNorm);... - FF(1,:,:);
            end
    end
end