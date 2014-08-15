InitBreach;
%% load formula
formulas = QMITL_ReadFile('spec.stl');
phi_template = phi_100;
  
%% Create system and input strategy 
mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'VarStep2');
Sys.tspan = 0:.01:50;
    

%% Property parameters 
prop_opt.params = {'t1', 'rpm_min'};
prop_opt.monotony = [-1 -1];
prop_opt.p_tol = [.1 5];
prop_opt.order = [1 2];
prop_opt.p_interval = [0 50; ...
                       0 6000; ... 
                       ];
  
%% System and falsification parameters
falsif_opt.params = {'dt_u0', 'throttle_u0', 'brake_u1'};
falsif_opt.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 325 ;
                       ];
  
falsif_opt.iter_max = 2;
falsif_opt.nb_init  = 2;
falsif_opt.nb_iter = 100;  


[p,rob] = ReqMining(Sys, phi_template, falsif_opt, prop_opt)
