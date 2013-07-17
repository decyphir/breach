
InitBreach;
%% load formula
formulas = QMITL_ReadFile('spec.stl');
phi_template = phi_stay2;

%% Create system and input strategy 
mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'VarStep2');
Sys.tspan = 0:.01:50;

%% Property parameters 
prop_opt = []; 
prop_opt.params = { 't1'};
prop_opt.monotony = [-1];
prop_opt.p_tol = [.02];
prop_opt.p_interval = [0 40]
               
%% System and falsification parameters
falsif_opt = [];
falsif_opt.params = {'dt_u0', 'throttle_u0', 'brake_u1'};
falsif_opt.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 325 ;
                       ];
  
falsif_opt.iter_max = 5;
falsif_opt.nb_init =  100;
falsif_opt.nb_iter = 0;  

[p,rob] = ReqMining(Sys, phi_template, falsif_opt, prop_opt)
