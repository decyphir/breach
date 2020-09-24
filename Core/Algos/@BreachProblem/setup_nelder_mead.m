function solver_opt = setup_nelder_mead(this,x0, varargin)


solver_opt = optimset();
solver_opt = optimset(solver_opt, 'Display', 'off');
this.x0 = (this.ub+this.lb)/2;

if nargin>=2
    this.x0 = x0;
end
if nargin>=3
    solver_opt = optimset(varargin{:});
end

this.solver = 'nelder_mead';
this.solver_options = solver_opt;


end
