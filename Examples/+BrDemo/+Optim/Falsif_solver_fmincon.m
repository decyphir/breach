%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
%% Choose the max number of cost function evaluations (each evaluation requires a simulation)
max_obj_eval = 500;

% Setup the solver using fmincon algorithm. Since fmincon is a 
% local search algorithm, we also apply multi-restart global search 
% methods to fminsearch. Thus, if fminsearch converges before the max
% cost evaluation allowed, it will restart the search from another random
% initial state. 
pb_se_fc = AFC_falsif03.copy();
pb_se_fc.max_obj_eval = max_obj_eval;
pb_se_fc.setup_solver('fmincon');
disp('============== fmincon with multiple restarts results ==============');
tic
[false_fc, res_fc] = pb_se_fc.solve();
toc