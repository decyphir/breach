%% Initialization
BrDemo.Optim.InitAFC_Falsif()
AFC_falsif03 = FalsificationProblem(B, phi03);
max_obj_eval = 500;
%% 
pb_se_ga = AFC_falsif03.copy();
pb_se_ga.max_obj_eval = max_obj_eval;
pb_se_ga.max_time = inf;
pb_se_ga.setup_solver('ga');
disp('============== genetic algorithm results ==============');
tic
[false_ga res_ga] = pb_se_ga.solve(); 
toc
%% 

pb_se_nm = AFC_falsif03.copy();
pb_se_nm.max_obj_eval = 100;
pb_se_nm.max_time = inf;
pb_se_nm.display = 'on';
disp('============== Nelder-Mead with multiple restarts results ==============');
tic
[false_nm res_nm] = pb_se_nm.solve();
toc

%%
pb_se_sa = AFC_falsif03.copy();
pb_se_sa.max_obj_eval = max_obj_eval;
pb_se_sa.max_time = inf;
pb_se_sa.setup_solver('simulannealbnd');
disp('============== simulannealbnd results ==============');
tic
[false_sa res_sa] = pb_se_sa.solve();
toc

%%
pb_se_fs = AFC_falsif03.copy();
pb_se_fs.max_obj_eval = max_obj_eval;
pb_se_fs.max_time = inf;
pb_se_fs.setup_solver('fminsearch');
disp('============== fminsearch with multiple restarts results ==============');
tic
[false_fs res_fs] = pb_se_fs.solve();
toc
%%
pb_se_fc = AFC_falsif03.copy();
pb_se_fc.max_obj_eval = max_obj_eval;
pb_se_fc.max_time = inf;
pb_se_fc.setup_solver('fmincon');
disp('============== fmincon with multiple restarts results ==============');
tic
[false_fc res_fc] = pb_se_fc.solve();
toc
%%
pb_se_cs = AFC_falsif03.copy();
pb_se_cs.max_obj_eval = max_obj_eval;
pb_se_cs.setup_solver('cmaes');
disp('============== cmaes results ==============');
tic
[false_cs res_cs] = pb_se_cs.solve();
toc

