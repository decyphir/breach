function res = solver_global_neldear_mead_legacy(this)

if this.solver_options.use_param_set_as_init
    X0 = this.BrSet.GetParam(this.params);
    this.BrSys.ResetParamSet;    % needs double-checking...
else
    X0 = init_basic_X0(this);
end

% display header
if ~strcmp(this.display,'off')   
    fprintf('Eval objective function on %d initial parameters.\n', size(X0,2));
    this.display_status_header();
end
res_init = FevalInit(this, X0);
res{1} = res_init;
res{1}.x0 = X0;
res{1}.fval = res_init.f;

this.solver_options.start_at_trial = this.solver_options.start_at_trial+this.solver_options.nb_new_trials;

if (this.solver_options.nb_local_iter>0) && (~this.stopping)
    fun_obj = @(x)(min(this.objective(x),[],1)); % for multi-objective support

    rfprintf_reset()
    if ~strcmp(this.display,'off')   
        fprintf('\nStarting local optimization using Nelder-Mead algorithm\n');
        this.display_status_header();
    end
    % Collect and sort solutions
    [~, ibest] = sort(max(res_init.fval,[],1));
    options = optimset(this.solver_options.local_optim_options, 'MaxIter',this.solver_options.nb_local_iter);
    flag_Cont = true;
    
    if this.use_parallel&&false % FIXME: behavior need be more thoroughly tested/verified 
        num_works = this.BrSys.Sys.Parallel;
        options = optimset(options, 'Display', 'off');
        options = optimset(options, 'UseParallel', true);
        while flag_Cont
            fun = @(x0) optimize(...
                fun_obj,x0,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
            for idx = 1:num_works
                x0 = X0(:,idx);
                F(idx) = parfeval(fun, 4, x0);
            end
            for idx = 1:num_works
                [completedIdx, x, fval, exitflag, output] = fetchNext(F); 
                res{end+1} = struct('x',x, 'fval',fval, 'exitflag', exitflag,  'output', output);
                % update the nb_obj_eval
                this.nb_obj_eval = this.nb_obj_eval + output.funcCount;
                this.X_log = [this.X_log output.logs];
                if fval < 0 || this.nb_obj_eval > this.max_obj_eval
                    cancel(F);
                    flag_Cont = false;
                    break;
                end
            end
            
        end 
    else

        for i_loc = ibest
            x0 = X0(:,i_loc);           
            if ~this.stopping()
                options = optimset(options, 'Display', 'off');
                [x, fval, exitflag, output] = optimize(...
                    fun_obj, x0 ,this.lb,this.ub,this.Aineq,this.bineq,this.Aeq,this.beq,[],[],options,'NelderMead');
                res{end+1} = struct('x0', x0, 'x',x, 'fval',fval, 'exitflag', exitflag,  'output', output);
            end
       end
    end
end

end
