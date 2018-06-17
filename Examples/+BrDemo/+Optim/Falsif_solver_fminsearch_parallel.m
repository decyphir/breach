%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Setup the solver using fminsearch algorithm
pb_par_fs = AFC_falsif03.copy();
pb_par_fs.max_obj_eval = max_obj_eval;
pb_par_fs.SetupParallel();
pb_par_fs.setup_solver('fminsearch');
% For fminsearch algorithm, we will run multiple optimization instances
% from different initial points.

% Due to the Mathworks' parallelism implemention, the log function provided 
% by BreachProblem won't be able to record the parameters tried during 
% optimization. However, the file cache function still works. 

% Since fminsearch is a local optimizer, breach wil try to bound the total 
% number of simulation allowed by max_obj_eval.

% The final result will return the best one among all parallel instances. 
disp('============== parallel fminsearch results ==============');
tic
[false_fs res_fs] = pb_par_fs.solve();
toc
pb_par_fs.StopParallel();