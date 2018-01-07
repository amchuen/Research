classdef boundCond
    
    properties
        % main roles
        vals % container for bc values... cell matrix?
        physTypes % container for physical bc types
        numTypes % container for numerical bc types
        ranges % container storing ranges for boundary conditions
        updateFlags % container for update flags (needed for dynamic BC's like patch's, non-flat walls, and symmetric BC's) 
        
        % grid data
        gr
        
        % control
        vecNum % corresponding fieldVec number
        numEqns % checks for number of equations solved (3 for isentropic, 4, for full Euler)
        isPolar % checks if calculations are in polar coords
    end
    
    methods % initialization routines
        % store bc's
        function obj = boundCond(caseStruct, gr, numEqns)
            % Note about caseStruct
            % assume all inputs will be put into cell matrix, even just
            % simple strings
            
            % Immediate initializations
            obj.gr = gr;
            obj.numEqns = numEqns;
            obj.vals = caseStruct.vals;
            obj.physTypes = caseStruct.physTypes;
            obj.ranges = caseStruct.ranges;
            
            % Manage mathtypes
            obj = obj.genNumTypes();
            
            % Check for updateFlags
            obj = obj.checkFlags(); 
            
        end
        
        function obj = genNumTypes(obj)
            DIR = fieldnames(obj.physTypes);
            for i = 1:length(DIR) % loop through each direction
                nBCs = length(obj.physTypes.(DIR{i})); % number of boundaries/direction
                obj.numTypes.(DIR{i}) = cell(obj.numEqns, nBCs);
                
                for ii = 1:nBCs                
                    switch obj.physTypes.(DIR{i}){ii}
                        case {'inlet', 'far-field', 'patch'}
                           obj.numTypes.(DIR{i}){:,ii} = 'D';
                        case 'outlet'
                           obj.numTypes.(DIR{i}){:,ii} = 'N';
                        case {'sym', 'wall'} % for the time being, it is assumed the wall behaves like a symmetry boundary, but grid may actually intersect w/ symmetry boundary
                            obj.numTypes.(DIR{i}){:,ii} = 'N';
                            if (strcmpi(DIR{i},'W') || strcmpi(DIR{i},'E'))
                                obj.numTypes.(DIR{i}){2,ii} = 'D';
                            elseif (strcmpi(DIR{i},'N') || strcmpi(DIR{i},'S'))
                                obj.numTypes.(DIR{i}){3,ii} = 'D';
                            end
                    end
                end
                
            end
            
        end
        
        function obj = checkFlags(obj)
            DIR = fieldnames(obj.physTypes);
            for i = 1:length(DIR) % loop through each direction
                nBCs = length(obj.physTypes.(DIR{i})); % number of boundaries/direction
                obj.updateFlags.(DIR{i}) = cell(obj.numEqns, nBCs);
                
                for ii = 1:nBCs                
                    switch obj.physTypes.(DIR{i}){ii}
                        case {'inlet', 'far-field', 'outlet', 'sym', 'wall'} % stays constant
                           obj.updateFlags.(DIR{i}){:,ii} = 0;
                        case 'patch' % needs to be updated with other Euler-Fields
                           obj.updateFlags.(DIR{i}){:,ii} = 1;
                    end
                end
                
            end
        end
        
    end
    
    methods 
        
        % update bc's for patches
        function obj = updateBC(obj, fieldObj, varargin)
            % Inputs:
            % obj -> calls self
            % FVobj -> used for getting boundary values if updateFlags
            % are triggered
            % varargin -> stores other objects that are needed?
            DIR = fieldnames(obj.updateFlags);
            
            for i = 1:length(DIR)
                nBCs = length(obj.physTypes.(DIR{i})); % number of boundaries/direction
                
                for ii = 1:nBCs 
                    flags = cell2mat(obj.updateFlags.(DIR{i}){:,ii});
                    
                    if any(flags == 1) % use the corresponding FVobj
                        if all(flags) % we are a vector normal to the wall
                           obj.vals 
                        else
                            
                            
                        end
                        
                    elseif any(flags == 2) % use 
                        
                        
                    end
                    
                end
                
            end
            
        end
        
        % return bc values
        function BCout = getBCval(obj, FVobj, DIR) 
            % OUTPUT: this function outputs the bc values for the given
            % direction... used for calculating bondary conditions
            nBCs = length(obj.physTypes.(DIR)); % number of boundaries/direction
            for i = 1:nBCs
                % extract the domain of the current boundary
                if (strcmpi(DIR,'W') || strcmpi(DIR,'E'))
                    domain = obj.gr.d2; 
                elseif (strcmpi(DIR,'N') || strcmpi(DIR,'S'))
                    domain = obj.gr.d1;
                end
                BCout = zeros([size(domain), obj.numEqns]);
                rangeInd = ((domain>obj.ranges.(DIR){i}(1)) & (domain<obj.ranges.(DIR){i}(2)));
                switch obj.physTypes.(DIR){i}
                    case {'Inlet', 'Outlet', 'Far-Field', 'Patch'}
                        BCout(rangeInd, :) = reshape(cell2mat(obj.vals.(DIR){:,i}),1,1,obj.numEqns);
                    case {'wall', 'sym'}
                        if (strcmpi(DIR,'W') || strcmpi(DIR,'E')) && (FVobj.vecNum == 2)
                            BCout(rangeInd) = % figure out how to get wall velocities here!
                        elseif (strcmpi(DIR,'N') || strcmpi(DIR,'S')) && (FVobj.vecNum == 3)
                            BCout(rangeInd) = 
                        end
                end
                
            end
            
        end
        
        
    end
            
    
end