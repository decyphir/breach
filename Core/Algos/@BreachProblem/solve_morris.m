function res = solve_morris(this)
% solve_corners works with quasi-random sampling
%

r = this.solver_options.num_path;
p = this.solver_options.size_grid;
dim = numel(this.params);

% number of samples is going to be (dim+1)*r
num_samples = (dim+1)*r;
if this.max_obj_eval<num_samples % we have to reevaluate the number of path
    r = floor(this.max_obj_eval/(dim + 1));
end

Sys= CreateSystem({},this.params, this.x0);
P  = CreateParamSet(Sys, this.params, [this.lb, this.ub]);
Pr = pRefine(P, p,r);
X0 = Pr.pts;

if ~strcmp(this.display,'off')
    fprintf('\n Running %g samples using Morris'' method\n', size(X0, 2));
    this.display_status_header();
end

res = this.FevalInit(X0);

[res.mu, res.mustar, res.sigma, res.sigmastar, res.EE] = EEffects(res.fval, Pr.D, p);
res.params = this.params;
this.add_res(res);

if ~strcmp(this.display, 'off')
    display_morris_result(res);
end

end

