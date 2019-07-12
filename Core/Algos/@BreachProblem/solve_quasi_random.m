function res = solve_quasi_random(this)
% solve_quasi_random works with quasi-random sampling   
% 

seed = this.solver_options.quasi_rand_seed;
num_samples =  this.solver_options.num_quasi_rand_samples;

if ~strcmp(this.display,'off')
    fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++\nRunning %g quasi-random samples with seed %g\n', num_samples, seed);
    this.display_status_header();
end

BrQ = this.BrSet.copy();
BrQ.ResetParamSet();
BrQ.SetParamRanges(this.params, [this.lb this.ub])
BrQ.QuasiRandomSample(num_samples, seed);
X0 = BrQ.GetParam(this.params);

res = this.FevalInit(X0);
this.add_res(res);
    

end
