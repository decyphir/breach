function res = solve_global_nelder_mead(this)

if this.solver_options.use_param_set_as_init
    X0 = this.BrSet.GetParam(this.params);
    this.BrSys.ResetParamSet;    % needs double-checking...
else
    X0 = init_basic_X0(this);
end

% display header
fprintf('Eval objective function on %d initial parameters.\n', size(X0,2));

this.display_status_header();
res = FevalInit(this, X0);
this.solver_options.start_at_trial = this.solver_options.start_at_trial+this.solver_options.nb_new_trials;

if (this.solver_options.nb_local_iter>0) && (~this.stopping)
    rfprintf_reset()
    fprintf('\nStarting local optimization using Nelder-Mead algorithm\n');
    
    this.display_status_header();
    % Collect and sort solutions
    [~, ibest] = sort(res.fval);
    options = optimset(this.solver_options.local_optim_options, 'MaxIter',this.solver_options.nb_local_iter);
    if this.search_parallel
        num_works = this.BrSys.Sys.Parallel;
        options = optimset(options, 'Display', 'off');
        %fun = @(x) 100*(x(1)/1000-0.95)^2 + (x(2)-20)^2 + (x(3)-37)^2;
        for idx = 1:num_works
            x0 = X0(:,idx);
            F(idx) = parfeval(@optimize, 4, ...
                this.objective, x0 ,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
            %fun, x0 ,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
        end
        res = cell(1, num_works);
        for idx = 1:num_works
            [completedIdx, sol, fval, exitflag, output] = fetchNext(F);
            res{completedIdx} = struct('sol',sol, 'fval',fval, 'exitflag', exitflag,  'output', output);
        end
    else
        for i_loc= ibest
            x0 = X0(:,i_loc);
            if ~this.stopping()
                [sol, fval, exitflag, output] = optimize(...
                    this.objective, x0 ,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
                res{end+1} = struct('sol',sol, 'fval',fval, 'exitflag', exitflag,  'output', output);
            end
        end
    end
end

end
