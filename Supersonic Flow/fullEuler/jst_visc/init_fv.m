function [BC, BC2] = init_fv(BC_type, BC_val, XX)

%% Solution

fnames_D = {'W', 'E', 'N', 'S'};
fnames_N = {'Wx', 'Ex', 'Ny', 'Sy'};

% BC_type is defined using W, E, N, S
% Different BC_types:
% * I -> inlet, dirichlet
% * O -> outlet, Neumann
% * W -> wall/reflecting, Neumann
% * S -> symmetry, Neumann

BC = struct();
BC2 = BC;

for i = 1:length(fnames_D)
    [r_check, c_check] = size(BC_val{i});
    if (r_check > 1) || (c_check > 1)
        % assume the value is user-defined
        tempBC = BC_val{i};
    else % single value BC
        if i >= 3 % check along x_dir for yBC
            tempBC = BC_val{i} * ones(1, size(XX, 2));
        else % check along y-dir for xBC
            tempBC = BC_val{i} * ones(size(XX, 1), 1);
        end
        
    end
    
    % Check for dirichlet or Neumann
    if strcmp(BC_type{i}, 'D') % Dirichlet condition occurs
        BC.(fnames_D{i}) = tempBC;
        BC2.(fnames_D{i}) = tempBC.^2;
    elseif strcmp(BC_type{i}, 'N') % Neumann condition occurs
        BC.(fnames_N{i}) = tempBC;
        BC2.(fnames_N{i}) = tempBC;
    end
    
end



end