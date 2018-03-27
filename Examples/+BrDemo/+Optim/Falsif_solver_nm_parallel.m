%% Initialization the falsification problems using AFC model
BrDemo.Optim.InitAFC_Falsif()
max_obj_eval = 500;

%% Use the default falsification algorithm Nelder-Mead
phi03 = set_params(AF_alw_ok, 'tol', 0.03);
pb_par_nm = FalsificationProblem(B, phi03);
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
pb_par_nm.StopParallel();