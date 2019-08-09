function solver_opt = setup_ga(this)
%  setup_ga genetic algorithm from global optimization toolbox
% 

this.solver= 'ga';
solver_opt = gaoptimset('TimeLimit',this.max_time, 'Display','off');

if numel(this.params) <= 5 
    pop_size = 50;
else
    pop_size = 200;
end

if isa(this, 'FalsificationProblem')
   if this.StopAtFalse
       solver_opt = gaoptimset(solver_opt, 'FitnessLimit',0);
   end
end

if isa(this, 'MaxSatProblem')
   if this.StopAtTrue
       solver_opt = gaoptimset(solver_opt, 'FitnessLimit',0);
   end
end

if this.use_parallel
    solver_opt = gaoptimset(solver_opt,'Vectorized', 'on'); 
end

% for multi-objective
% https://www.mathworks.com/help/gads/gamultiobj.html

% for 2018 version matlab using OPTIMOPTIONS and 
% set 'MaxGenerations' to 100 * length(this.params));

solver_opt = gaoptimset(solver_opt, 'Generations', 100 * length(this.params));

this.solver_options = solver_opt;

end
