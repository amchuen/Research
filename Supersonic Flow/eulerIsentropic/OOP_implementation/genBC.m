function BC = genBC(phys_types, vals, range, varargin)
% inputs:
    % - phys_types: physical boundaries of grid, each direction set should be
    % enclosed in cell for each element of input cell
    % - vals: physical values corresponding to dirchlet/neumann assignments... defaults are assigned if empty
    % - range: sets range for each type present at each direction... empty
    % implies that the direction only has one type

BC = struct('W', struct, 'N', struct, 'E', struct, 'S', struct);
for vec = 1:length(vals{1}{1})
   BC(vec) = BC(1); 
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
        BC(vec).(dir{i}).physical = phys_types{i};
        if (isempty(range)|| isempty(range{i})) && (n_types == 1) % assumes that there is only one BC at each direction
            BC(vec).(dir{i}).range = {};
        elseif (n_types > 1)
            BC(vec).(dir{i}).range = range{i};
        else
            error('Please enter ranges!\n');
        end
        
        for ii = 1:n_types % loop through types at each direction
            % default properties
            BC(vec).(dir{i}).update{ii} = 0;

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
                    BC(vec).(dir{i}).update{ii} = ~all(vals{i}{ii}{wall_vec} == 0);

                    if vec == wall_vec
                        BC(vec).(dir{i}).numerical{ii} = repmat({'D'}, 1, length(vals{1}{1}));
                        BC(vec).(dir{i}).val{ii} = repmat({vals{i}{ii}{vec}}, 1, length(vals{1}{1}));
                        BC(vec).(dir{i}).val{ii}{vec} = BC(vec).(dir{i}).val{ii}{vec}.*vals{i}{ii}{vec};

                    else
                        BC(vec).(dir{i}).numerical{ii} = repmat({'N'}, 1, length(vals{1}{1}));
                        BC(vec).(dir{i}).numerical{ii}{wall_vec} = 'D';
                        BC(vec).(dir{i}).val{ii} = vals{i}{ii};
                    end

                case 'inlet'
                    BC(vec).(dir{i}).numerical{ii} = repmat({'D'}, 1, length(vals{1}{1}));
                    BC(vec).(dir{i}).val{ii} = vals{i}{ii};

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       BC(vec).(dir{i}).val{ii}{iii} = BC(vec).(dir{i}).val{ii}{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
                    end

                case 'outlet'
                    BC(vec).(dir{i}).numerical{ii} = repmat({'N'}, 1, length(vals{1}{1}));
                    BC(vec).(dir{i}).val{ii} = {0,0,0};

                case 'far-field'
                    BC(vec).(dir{i}).numerical{ii} = repmat({'D'}, 1, length(vals{1}{1}));
                    BC(vec).(dir{i}).val{ii} = vals{i}{ii};

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       BC(vec).(dir{i}).val{ii}{iii} = BC(vec).(dir{i}).val{ii}{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
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
                    
                    BC(vec).(dir{i}).numerical{ii} = repmat({'N'}, 1, length(vals{1}{1}));
                    BC(vec).(dir{i}).numerical{ii}{sym_vec} = 'D';

                case 'patch'
                    BC(vec).(dir{i}).update{ii} = 1;
                    BC(vec).(dir{i}).numerical{ii} = repmat({'D'}, 1, length(vals{1}{1}));
                    BC(vec).(dir{i}).val{ii} = vals{i}{ii}; % initial values given assumed to be similar to initial condition values

                    for iii = 1:length(vals{i}{ii}) % loop through elements of vector
                       BC(vec).(dir{i}).val{ii}{iii} = BC(vec).(dir{i}).val{ii}{iii} .* vals{i}{ii}{vec} ./ vals{i}{ii}{1};
                    end
            end
        end

    end

end


end