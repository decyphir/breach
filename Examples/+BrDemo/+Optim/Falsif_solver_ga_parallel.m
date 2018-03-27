%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Setup the solver using genetic algorithm
pb_par_ga = AFC_falsif03.copy();
pb_par_ga.max_obj_eval = max_obj_eval;
pb_par_ga.SetupParallel();
pb_par_ga.setup_solver('ga');
% For genetic algorithm, since it is a global optimizer, breach chooses to
% parallelize the children exploration but with one optimization instance.
disp('============== parallel genetic algorithm results ==============');
tic
[false_ga res_ga] = pb_par_ga.solve(); 
toc
pb_par_ga.StopParallel();