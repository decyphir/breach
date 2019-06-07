function res = solve_binsearch(this)
% solve_binsearch 

Sys = this.BrSys.Sys;
phi = STL_Formula(this.Spec.req_monitors{1}.formula_id);
P = this.BrSet.P;
opt = this.solver_options;

% Infer monotonicity
if strcmp(opt.monotony, 'infer')
    Br = this.BrSet.copy();
    if this.solver_options.verbose
        fprintf('Inferring monotonicity...\n')
    end
    opt.monotony = Br.ChecksMonotony(this.Spec, this.params, [this.lb this.ub]);
    if any(opt.monotony==0)
        error('solve_binsearch:not_monotonic','The problem does not appear to be monotonic wrt parameters, consider setting it manually in this.solver_options.monotony (+1 if inc, -1 if dec) or try another solver.');
    else
        if this.solver_options.verbose
            fprintf('The problem appears to be monotonic in\n')
            for ip=1:numel(this.params)
                if opt.monotony(ip)==1
                    fprintf([this.params{ip} ' (increasing)\n'])
                else
                    fprintf([this.params{ip} ' (decreasing)\n'])
                end
            end
            fprintf('\n');
        end
        this.solver_options.monotony = opt.monotony;
    end
end

[p, rob] = GetPropParamBin(Sys, phi, P, opt,P.traj, this.T_Spec);
res.x = p;
res.f = rob;
this.x_best = p;
this.obj_best = rob;
this.X_log = p;
this.obj_log = rob;
this.BrSet.SetParam(this.params, p, true);
this.BrSet.CheckSpec(phi);
end
