function res = solve_ga(this, problem)

% init seed if needed
if isfield(this.solver_options,'rand_seed')
    rng(this.solver_options.rand_seed);
end

this.display_status_header();
                                
[x,fval,exitflag,output, population, scores] = ga(problem);
res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);

% TODO ? provide LogX as Output function for ga, run after each generation
%if this.use_parallel
%    this.LogX(population', scores');
%end

this.add_res(res);

end
