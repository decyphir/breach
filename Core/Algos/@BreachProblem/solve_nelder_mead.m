function res = solve_nelder_mead(this)

x0 = this.x0;
    
% display header
if ~strcmp(this.display,'off')
    fprintf('\n\n********************************************\nStarting Nelder Mead optimization from x0:\n');
    this.Display_X(x0);
    this.display_status_header();
end

opt = this.solver_options;

fun_obj = @(x)(min(this.objective(x),[],1)); % for multi-objective support
[x, fval, exitflag, output] = minimize(...
    fun_obj, x0 ,this.Aineq, this.bineq,this.Aeq, this.beq, this.lb,this.ub,[],opt);

res = struct('x0', x0, 'x',x, 'fval',fval, 'exitflag', exitflag,  'output', output);
this.add_res(res);

end

