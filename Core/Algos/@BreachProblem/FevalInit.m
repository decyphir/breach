function res = FevalInit(this,X0)
% FevalInit Eval objective function on parameters in BrSet

if ~exist('X0', 'var')
    X0 = this.BrSet.GetParam(this.params);
end
fval = this.objective(X0); 
res = struct('x',this.x_best,'f', this.obj_best, 'fval', fval);
this.res = res;

end

