%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Setup the solver using cmaes algorithm
pb_par_cs = AFC_falsif03.copy();
pb_par_cs.max_obj_eval = max_obj_eval;
% For the cmaes solver, the parallelism is implemented by running children 
% exploraion in parallel which means falsificaiton has only one search 
% instance but simulations needed during children exploraion run in parallel)
pb_par_cs.SetupParallel();
pb_par_cs.setup_solver('cmaes');
disp('============== parallel cmaes results ==============');
tic
[false_cs res_cs] = pb_par_cs.solve();
toc
pb_par_cs.StopParallel();