function res = solve_ga(this, problem)

[x,fval,exitflag,output, population, scores] = ga(problem);
res =struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);

this.x_best = x;
this.obj_best = fval;

if this.use_parallel
    this.LogX(population', scores);
end

this.res=res;

end
