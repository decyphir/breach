function monotony = check_monotony(this)
% BreachProblem.check_monotony might evolve eventually into more general sensitivity
% analysis

% check domains for non double -
for id= 1:numel(this.domains)
   if ~isequal(this.domains(id).type, 'double')
       warning('check_monotony not tested with non double type, such as type %s of variable %s', this.domains(id), this.params{id});
   end
end

Sys= CreateSystem({},this.params, this.x0);
P  = CreateParamSet(Sys, this.params, [this.lb, this.ub]);

p = 4; 
r = 10; 
Pr = pRefine(P, p,r);
Y = this.objective(Pr.pts);
[~, ~, ~, ~, EE] = EEffects(Y, Pr.D, p);

monotony = all(EE'>=0)-all(EE'<=0); % 1 if all positive, -1 if all negative, 0 otherwise



