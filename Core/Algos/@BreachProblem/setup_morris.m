function solver_opt = setup_morris(this, varargin)
% solve_corners works with quasi-random sampling   
% 

solver_opt = struct('num_path', 100, ...   % number of paths, i.e., set of
                    ...                    % samples providing 1 pair of samples per dim
                'rand_seed', 1, ...         % random seed for reproducibility    
                'size_grid', 5 ...     % number of grid levels, intervals in each dim
    );
solver_opt = varargin2struct_breach(solver_opt, varargin{:});


this.solver = 'morris';
this.solver_options = solver_opt;                      

end
