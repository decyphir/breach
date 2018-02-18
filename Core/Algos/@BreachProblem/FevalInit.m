function res = FevalInit(this,X0)
% FevalInit Eval objective function on parameters in BrSet

if ~exist('X0', 'var')
    X0 = this.BrSet.GetParam(this.params);
end

fval = this.objective(X0); 

[fopt,iopt] = min(fval);
xopt = X0(:,iopt);

if fopt < this.obj_best
    this.x_best = xopt;
    this.obj_best = fopt;
end

res = struct('x',xopt,'f', fopt, 'fval', fval);
this.res = res;

end

