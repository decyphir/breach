%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
%% Choose the max number of cost function evaluations (each evaluation requires a simulation)
max_obj_eval = 500;

%% Setup the solver using Nelder-Mead algorithm
pb_se_nm = AFC_falsif03.copy();
pb_se_nm.max_obj_eval = max_obj_eval;
pb_se_nm.display = 'on';
% the default solver in Breach is Nelder-Mead (no need to setup the solver)
disp('============== Nelder-Mead with multiple restarts results ==============');
tic
[false_nm res_nm] = pb_se_nm.solve();
toc