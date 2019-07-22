function res = FevalInit(this,X0)
% FevalInit Eval objective function on parameters in BrSet

if ~exist('X0', 'var')
    X0 = this.BrSet.GetParam(this.params);
end
[fval, cval] = this.objective(X0);

[fbest, ibest] = min(min(fval));
res = struct('BrSys', this.BrSys, 'X0',X0,'x',X0(:,ibest),'f', fbest, 'fval', fval,'cval',cval);

end

