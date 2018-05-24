%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Setup the solver using fmincon algorithm
pb_par_fc = AFC_falsif03.copy();
pb_par_fc.max_obj_eval = max_obj_eval;
pb_par_fc.SetupParallel();
pb_par_fc.setup_solver('fmincon');
% For fmincon?algorithm, we use the Mathworks' native parallelism support 
% 'UseParallel' flag for the optimoptions setting
disp('============== parallel fmincon with multiple restart results ==============');
tic
[false_fc res_fc] = pb_par_fc.solve();
toc
pb_par_fc.StopParallel();