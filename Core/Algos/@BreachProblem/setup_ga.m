function solver_opt = setup_ga(this)
%  setup_ga genetic algorithm from global optimization toolbox
% 

this.solver= 'ga';
this.display = 'off'; % use ga display by default

%solver_options = gaoptimset('TimeLimit',this.max_time);
solver_options = gaoptimset('TimeLimit',this.max_time,... % doesn't seem to work well in parallel
    'Display', 'iter');


if this.max_obj_eval<inf
    solver_options = gaoptimset('Generations', ceil(this.max_obj_eval/50)); % check that 50 is default population size..
end

if isa(this, 'FalsificationProblem')
   if this.StopAtFalse
    solver_options = gaoptimset(solver_options, 'FitnessLimit',0); 
   end
end
if isa(this, 'MaxSatProblem')
   if this.StopAtTrue
    solver_options = gaoptimset(solver_options, 'FitnessLimit',0); 
   end
end

if this.use_parallel
    solver_options= gaoptimset(solver_options,'UseParallel', true);
end

this.solver_options = solver_options;

end
