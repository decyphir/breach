function opt = setup_global_nelder_mead(this, varargin)

this.solver = 'global_nelder_mead';
dim_pb = numel(this.params);
opt = struct( ...
    'num_corners', min(2^dim_pb, 10*dim_pb),...
    'group_by_inputs', true,...
    'num_quasi_rand_samples', 10*dim_pb, ...
    'quasi_rand_seed', 1, ...
    'local_max_obj_eval', this.max_obj_eval/2, ...
    'local_solver', 'nelder_mead',... 
    'local_options', optimset('Display', 'off') ...
    );


if this.use_parallel
    opt.use_parallel = true;
end

% checks with what we have already
if isstruct(this.solver_options)
    fn = fieldnames(this.solver_options);
    for ifn = 1:numel(fn)
        field = fn{ifn};
        if isfield(opt, field)
            opt.(field) = this.solver_options.(field);
        end
    end
end

if nargin>2    
    opt = varargin2struct_breach(opt, varargin{:});

end

this.solver_options = opt;

end
