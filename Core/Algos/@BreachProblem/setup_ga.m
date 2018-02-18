function solver_opt = setup_ga(this)
%  setup_ga genetic algorithm from global optimization toolbox
% 

this.solver= 'ga';
this.display = 'off'; % use ga display by default

solver_options = gaoptimset('TimeLimit',this.max_time,... 
    'Display', 'iter');

if numel(this.params) <=5 
    pop_size = 50;
else
    pop_size = 200;
end

if this.max_obj_eval<inf
    solver_options = gaoptimset('Generations', ceil(this.max_obj_eval/pop_size)); % check that 50 is default population size..
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

solver_options= gaoptimset(solver_options,'Vectorized', 'on');  % makes it work in parallel
this.solver_options = solver_options;

end
