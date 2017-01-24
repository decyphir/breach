function res = FevalInit(this,X0)
% FEvaliInit Eval objective function on parameters in BrSys

if ~exist('X0', 'var')
    X0 = this.BrSet.GetParam(this.params);
end
nb_init = min(size(X0,2), this.max_obj_eval);
fval = inf(1,nb_init);

if this.use_parallel
    
    fun = @(x) this.objective_fn(x);
    
    % Display header
    this.display_status_header();
    fq = this.freq_update;

    for iter=1:floor(nb_init/fq)
        
        i_range =(iter-1)*fq+1:fq*iter;
        parfor isample = i_range
            fval(isample) = fun(X0(:,isample));
        end
        this.LogX(X0(:, i_range), fval(i_range));
        
        % update status
        this.display_status();

        if this.stopping()
            break
        end
    end
    
    if ~this.stopping()
        i_range = fq*iter:nb_init;
        parfor isample = i_range
            fval(isample) = fun(X0(:,isample));
        end
        this.LogX(X0(:, i_range), fval(i_range));
        this.display_status();
    end
    
else % serial  
    fun = this.objective;
    for isample = 1:nb_init
        fval(isample) = fun(X0(:,isample));
    end
end

[fopt,iopt] = min(fval);
xopt = X0(:,iopt);

this.x_best = xopt;
this.obj_best = fopt;

%           res = struct('x',this.x_best,'f', this.obj_best,'fval', this.obj_log);
res = struct('x',xopt,'f', fopt, 'fval', fval);
this.res = res;

end

