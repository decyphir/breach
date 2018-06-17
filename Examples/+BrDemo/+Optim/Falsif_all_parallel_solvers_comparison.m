%% Initialization
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;
%% 
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
%%
% parallelism of the ga solver is using vectorization, the logX works 
size(pb_par_ga.X_log)
%%
pb_par_fc = AFC_falsif03.copy();
pb_par_fc.max_obj_eval = max_obj_eval;
pb_par_fc.SetupParallel();
pb_par_fc.setup_solver('fmincon');
% For fmincon algorithm, we use the Mathworks' native parallelism support 
% 'UseParallel' flag for the optimoptions setting
disp('============== parallel fmincon with multiple restart results ==============');
tic
[false_fc res_fc] = pb_par_fc.solve();

%%
% parallelism of the ga solver is using native parallel support from Mathowrks
% the logX works 
size(pb_par_fc.X_log)


%%
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
%%
% The logX function does not work
size(pb_par_sa.X_log)
%%
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

%%
% The logX function does not work
size(pb_par_fs.X_log)
%% 
phi03 = set_params(AF_alw_ok, 'tol', 0.03);
pb_par_nm = FalsificationProblem(B, phi03);
pb_par_nm.SetupDiskCaching();
pb_par_nm.max_obj_eval = 100;
pb_par_nm.display = 'off';
pb_par_nm.SetupParallel();
% For Nelder-Mead, we will run multiple optimization instances from
% promising initial points (after the initial sample phase)

% Breach will try to respect the bound on the total number of simulation
% allowed by max_obj_eval.
disp('============== parallel Nelder-Mead with multiple restart results ==============');
tic
[false_nm res_nm] = pb_par_nm.solve();
toc
%%
% The function returning the tried parameters is implemented by changing 
% the interface of NM solver engine
% However, the number of traces in the cached folder is way more than the
% one logged.
size(pb_par_nm.X_log)
%%
pb_par_cs = AFC_falsif03.copy();
pb_par_cs.max_obj_eval = max_obj_eval;
pb_par_cs.SetupDiskCaching();
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
%%
% The logX function work however, the number of caches traces are more than
% logged. 
size(pb_par_cs.X_log)