%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
%% Choose the max number of cost function evaluations (each evaluation requires a simulation)
max_obj_eval = 500;

%% Setup the solver using genetic algorithm
pb_se_cs = AFC_falsif03.copy();
pb_se_cs.max_obj_eval = max_obj_eval;
pb_se_cs.setup_solver('cmaes');
disp('============== cmaes results ==============');
tic
[false_cs res_cs] = pb_se_cs.solve();
toc