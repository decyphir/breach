function solver_opt = setup_binsearch(this)
% setup_binsearch 
% TODO documentation

this.solver = 'binsearch';
this.display = 'off';
solver_opt = struct('params', {this.params},...
    'monotony', 'infer',...
    'ranges', [this.lb this.ub],...
    'verbose', 1);
this.solver_options = solver_opt;
end
