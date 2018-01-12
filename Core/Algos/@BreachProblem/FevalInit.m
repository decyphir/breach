function res = FevalInit(this,X0)
% FEvalInit Eval objective function on parameters in BrSys

if ~exist('X0', 'var')
    X0 = this.BrSet.GetParam(this.params);
end
nb_init = min(size(X0,2), this.max_obj_eval);
fval = inf(1,nb_init);

if this.use_parallel==0
    fun = this.objective;
    for isample = 1:nb_init
        fval(isample) = fun(X0(:,isample));
    end
else 
    fun = @(x) this.objective_fn(x);
    
    % Launch tasks 
    for idx = 1:nb_init
        par_f(idx) = parfeval(@(isample) fun(X0(:,isample)),1, idx);
    end
   
    % Display header
    this.display_status_header();
    fq = this.freq_update;

    for iter=1:nb_init
        [idx, value] = fetchNext(par_f);
        fval(idx) = value;
        this.LogX(X0(:, idx), fval(idx));
        
        % update status
        if rem(iter,fq)==0
            this.display_status();
        end
        if this.stopping()
            break
        end
    end
    cancel(par_f);
end

[fopt,iopt] = min(fval);
xopt = X0(:,iopt);

if fopt < this.obj_best
    this.x_best = xopt;
    this.obj_best = fopt;
end

%           res = struct('x',this.x_best,'f', this.obj_best,'fval', this.obj_log);
res = struct('x',xopt,'f', fopt, 'fval', fval);
this.res = res;

end

