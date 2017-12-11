function res = solve_global_nelder_mead(this)

if this.solver_options.use_param_set_as_init
    X0 = this.BrSet.GetParam(this.params);
    this.BrSys.ResetParamSet;    % needs double-checking...
else
    X0 = init_basic_X0(this);
end

% display header
fprintf('Eval objective function on %d initial parameters.\n', size(X0,2));

res = FevalInit(this, X0);
this.solver_options.start_at_trial = this.solver_options.start_at_trial+this.solver_options.nb_new_trials;

if (this.solver_options.nb_local_iter>0) && (~this.stopping)
    rfprintf_reset()
    fprintf('Local optimization using Nelder-Mead algorithm\n');
    
    % Collect and sort solutions
    [~, ibest] = sort(res.fval);
    options = optimset(this.solver_options.local_optim_options, 'MaxIter',this.solver_options.nb_local_iter);
    for i_loc= ibest
        x0 = X0(:,i_loc);
        if ~this.stopping()
            optimize(@(x)this.objective_wrapper(x),x0 ,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
        end
    end
end

end
