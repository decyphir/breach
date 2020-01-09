function res = solve_random(this)
% Solve_quasi_random works with quasi-random sampling   
% 

seed = this.solver_options.rand_seed;
num_samples =  this.solver_options.num_rand_samples;

if ~strcmp(this.display,'off')
    fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++\nRunning %g random samples with seed %g\n', num_samples, seed);
    this.display_status_header();
end

BrR = this.BrSet.copy();
BrR.ResetParamSet();
BrR.SetParamRanges(this.params, [this.lb this.ub])
rng(seed);
BrR.SampleDomain(this.params,num_samples);
X0 = BrR.GetParam(this.params);

res = this.FevalInit(X0);
this.add_res(res);
    
end
