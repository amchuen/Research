function varargout = manVars(task, varargin)

    if strcmpi(task, 'init') 
        % Assume we are creating new field variables
        % order of varargin: {ranges, dx, gam, M0, isPolar}
        ranges = varargin{1};
        dx = varargin{2};
        gam = varargin{3};
        M0 = varargin{4};
        isPolar = varargin{5};
        
        GR = struct('range1',[], 'range2',[], 'D1',[], 'D2',[]);
        EE = struct('fv',[],'f1',[], 'f2',[], 'f1_f',[], 'f1_b',[], 'f2_f',[], 'f2_b',[]);
        FF = EE; 
        GG = EE; 

        GR.range1 = ranges(2,1):dx(2):ranges(2,2);
        GR.range2 = ranges(1,1):dx(1):ranges(1,2);
        [GR.D2, GR.D1] = meshgrid(GR.range2, GR.range1);
        if isPolar
            GR.XX = GR.D1 .* cos(GR.D2);
            GR.YY = GR.D1 .* sin(GR.D2);
        else
            GR.XX = GR.D2;
            GR.YY = GR.D1;
        end

        % Matrix Dimensions
        % 1:Y, 2:X, 3:vec, 4:t

        EE.fv = repmat(cat(3,   ones(size(GR.D1)),... % rho
                                ones(size(GR.D1)),... % rho * u
                                zeros(size(GR.D1)),...% rho * v
                                ones(size(GR.D1))./(gam*(gam-1)*M0^2) + 0.5,... % rho * E
                                1/(gam * M0^2).*ones(size(GR.D1))),... % P ... eps_mach error here?
                                1, 1, 1, 3); 

        % EE.fv(:,:,5,:) = repmat((EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1),1,1,1,3);
        EE.fx_f = zeros(size(EE.fv(:,:,:,2)));
        EE.fx_b = EE.fx_f;
        EE.fy_f = EE.fx_f;
        EE.fy_b = EE.fx_f;

        FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
        FF.fx = zeros(size(FF.fv));

        GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
        GG.fy = zeros(size(GG.fv));
        
        % Assign Outputs
        varargout{1} = EE;
        varargout{2} = FF;
        varargout{3} = GG;
        varargout{4} = GR;
    
    elseif strcmpi(task, 'update')
        % Reassign inputs accordingly
        EE = varargin{1};
        FF = varargin{2};
        GG = varargin{3};
        gam = varargin{4};
        
        EE.fv(:,:,:,1:2) = EE.fv(:,:,:,2:3);
        EE.fv(:,:,5,2) = (EE.fv(:,:,4,2) - 0.5.*(EE.fv(:,:,2,2).^2 + EE.fv(:,:,3,2).^2)./EE.fv(:,:,1,2)).*(gam - 1);
        FF.fv = repmat(EE.fv(:,:,2,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
        GG.fv = repmat(EE.fv(:,:,3,2)./EE.fv(:,:,1,2),1,1,size(EE.fv,3)) .* EE.fv(:,:,:,2);
        
        % Assign Outputs
        varargout{1} = EE;
        varargout{2} = FF;
        varargout{3} = GG;
    end

end