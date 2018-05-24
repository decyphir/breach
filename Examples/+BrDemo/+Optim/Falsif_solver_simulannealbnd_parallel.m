%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Setup the solver using simulannealbnd algorithm
pb_par_sa = AFC_falsif03.copy();
pb_par_sa.max_obj_eval = max_obj_eval;
pb_par_sa.SetupParallel();
pb_par_sa.setup_solver('simulannealbnd');
% For simulannealbnd algorithm, we will run multiple global optimization
% instances in parallel. 

% Due to the Mathworks' parallelism implementation,
% the log function provided by BreachProblem won't be able to record the
% parameters tried during optimization. However, the file cache function 
% still works. 

% For each optimization instance, the number of simulation allowed is
% bounded by max_obj_eval.

% The final result will return the best one among all parallel instances.  
disp('============== parallel simulannealbnd results ==============');
tic
[false_sa res_sa] = pb_par_sa.solve();
toc
pb_par_sa.StopParallel();