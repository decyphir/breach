function solver_opt = setup_basic(this)
%  setup_basic basic solver: works with corners + quasi-random sampling   
% 

this.solver= 'basic';
solver_opt = struct( ...
    'start_at_trial', 0, ...
    'nb_new_trials', 10*numel(this.params) ...
    );
this.solver_options = solver_opt;
end
