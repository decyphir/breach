%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
%% Choose the max number of cost function evaluations (each evaluation requires a simulation)
max_obj_eval = 500;

%% Setup the solver using simulannealbnd algorithm
pb_se_sa = AFC_falsif03.copy();
pb_se_sa.max_obj_eval = max_obj_eval;
pb_se_sa.setup_solver('simulannealbnd');
disp('============== simulannealbnd results ==============');
tic
[false_sa res_sa] = pb_se_sa.solve();
toc