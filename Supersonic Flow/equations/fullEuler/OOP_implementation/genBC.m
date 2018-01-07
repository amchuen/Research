    function BC = genBC(phys_types, vals, range, varargin)
% inputs:
    % - phys_types: physical boundaries of grid, each direction set should be
    % enclosed in cell for each element of input cell
    % - vals: physical values corresponding to dirchlet/neumann assignments... defaults are assigned if empty
    % - range: sets range for each type present at each direction... empty
    % implies that the direction only has one type

BC = struct('W', struct, 'N', struct, 'E', struct, 'S', struct);
for vec = 1:length(vals{1}{1})
   obj.BC(vec) = obj.BC(1); 
end
% each dir should have the following fields
    % - physical type
    % - numerical type
    % - update ctrl
    % - values
    % - ranges

dir = fieldnames(BC);

for i = 1:length(dir) % loop through directions
    n_types = length(phys_types{i});

    for vec = 1:length(vals{1}{1}) % loop through the vectors
        
        for ii = 1:n_types % loop through types at each direction
            % default properties
            obj.BC(vec).(dir{i})(ii).update{ii} = 0;
            obj.BC(vec).(dir{i})(ii).physical = phys_types{i}{ii};
            if (isempty(range)|| isempty(range{i})) && (n_types == 1) % assumes that there is only one BC at each direction
                obj.BC(vec).(dir{i})(ii).range = {};
            elseif (n_types > 1)
                obj.BC(vec).(dir{i})(ii).range = range{i}{ii};
            else
                error('Please enter ranges!\n');
            end

            switch phys_types{i}{ii}
                case 'wall'          
                    % determine which velocity is normal to the
                    % wall
                    if (strcmp(dir{i},'W') || strcmp(dir{i}, 'E'))
                        wall_vec = 2;
                    elseif (strcmp(dir{i},'S') || strcmp(dir{i}, 'N'))
                        wall_vec = 3;
                    end

                    % check if wall is zero or a linearized BC
                    obj.BC(vec).(dir{i})(ii).update{ii} = ~all(vals{i}{ii}{wall_vec} == 0);

                    if vec == wall_vec
                        obj.BC(vec).(dir{i})(ii).numerical = repmat({'D'}, 1, length(vals{1}{1}));
                        obj.BC(vec).(dir{i})(ii).val= repmat({vals{i}{ii}{vec}}, 1, length(vals{1}{1}));
                        obj.BC(vec).(dir{i})(ii).val{vec} = obj.BC(vec).(dir{i})(ii).val{vec}.*vals{i}{ii}{vec};

                    else
                        obj.BC(vec).(dir{i})(ii).numerical = repmat({'N'}, 1, length(vals{1}{1}));
                        obj.BC(vec).(dir{i})(ii).numerical{wall_vec} = 'D';
                        obj.BC(vec).(dir{i})(ii).val = vals{i}{ii};
                    end

                case 'inlet'
                    obj.BC(vec).(dir{i})(ii).numerical = repmat({'D'}, 1, length(vals{1}{1}));
                    obj.BC(vec).(dir{i})(ii).val = vals{i}{ii};

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       obj.BC(vec).(dir{i})(ii).val{iii} = obj.BC(vec).(dir{i})(ii).val{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
                    end

                case 'outlet'
                    obj.BC(vec).(dir{i})(ii).numerical = repmat({'N'}, 1, length(vals{1}{1}));
                    obj.BC(vec).(dir{i})(ii).val = {0,0,0};

                case 'far-field'
                    obj.BC(vec).(dir{i})(ii).numerical = repmat({'D'}, 1, length(vals{1}{1}));
                    obj.BC(vec).(dir{i})(ii).val = vals{i}{ii};

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       obj.BC(vec).(dir{i})(ii).val{iii} = obj.BC(vec).(dir{i})(ii).val{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
                    end
                    
                case 'sym'
                    % in this case, if the numerical type is
                    % Neumann, then it is implied that the value is
                    % not normal to the surface, wherewas Dirichlet
                    % will be; i.e. the two distinctions will serve
                    % a different purpose

                    % need to check if boundary is on axis or not?
                    % updates not necessary, as can be extrapolated
                    % from the field values
                    if (strcmp(dir{i},'W') || strcmp(dir{i}, 'E'))
                        sym_vec = 2;
                    elseif (strcmp(dir{i},'S') || strcmp(dir{i}, 'N'))
                        sym_vec = 3;
                    end
                    
                    obj.BC(vec).(dir{i})(ii).numerical = repmat({'N'}, 1, length(vals{1}{1}));
                    obj.BC(vec).(dir{i})(ii).numerical{sym_vec} = 'D';

                case 'patch'
                    obj.BC(vec).(dir{i})(ii).update{ii} = 1;
                    obj.BC(vec).(dir{i})(ii).numerical = repmat({'D'}, 1, length(vals{1}{1}));
                    obj.BC(vec).(dir{i})(ii).val = vals{i}{ii}; % initial values given assumed to be similar to initial condition values

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       obj.BC(vec).(dir{i})(ii).val{iii} = obj.BC(vec).(dir{i})(ii).val{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
                    end
            end
        end

    end

end


end
