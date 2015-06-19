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
%prop_opt.params = {'rpm_min', 't1'};
prop_opt.monotony = [-1 -1];
prop_opt.p_tol = [.1 5];
prop_opt.order = [2 1];
prop_opt.ranges = [0 50; ...
                       0 6000; ... 
                       ];
  
%% System and falsification parameters
falsif_opt.params = {'dt_u0', 'throttle_u0', 'brake_u1'};
falsif_opt.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 325 ;
                       ];
  
falsif_opt.nb_init  = 2;
falsif_opt.nb_iter = 50;  
falsif_opt.nb_max_call = 100;

%% Max number of mining iterations
iter_max= 100;
[p,rob,Pr_speed_100_rpm] = ReqMining(Sys, phi_template, falsif_opt, prop_opt, iter_max);

Psave(Sys, 'Pr_speed_100_rpm');, % run Breach(Sys) to explore the successive counter-examples 
