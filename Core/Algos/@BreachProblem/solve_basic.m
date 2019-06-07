function res = solve_basic(this)
%  solve_basic basic solver: works with corners + quasi-random sampling   
% 
X0 = this.init_basic_X0();
this.display_status_header();
res = this.FevalInit(X0);
this.add_res(res); 
this.solver_options.start_at_trial = this.solver_options.start_at_trial+this.solver_options.nb_new_trials;

end
