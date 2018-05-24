%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
%% Choose the max number of cost function evaluations (each evaluation requires a simulation)
max_obj_eval = 500;

%% Setup the solver using genetic algorithm
pb_se_ga = AFC_falsif03.copy();
pb_se_ga.max_obj_eval = max_obj_eval;
pb_se_ga.setup_solver('ga');
disp('============== genetic algorithm results ==============');
tic
[false_ga res_ga] = pb_se_ga.solve(); 
toc