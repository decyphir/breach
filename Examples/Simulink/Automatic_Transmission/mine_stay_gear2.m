
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
prop_opt.ranges = [0 40]
               
%% System and falsification parameters
falsif_opt = [];
falsif_opt.params = {'dt_u0', 'throttle_u0', 'brake_u1'};
falsif_opt.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 325 ;
                       ];
 
falsif_opt.nb_init =  1000;
falsif_opt.nb_iter = 1000;  
falsif_opt.nb_max_call = 10000;

iter_max= 100;
[p,rob,Pr_stay_gear2] = ReqMining(Sys, phi_template, falsif_opt, prop_opt,iter_max)
Psave(Sys, 'Pr_stay_gear2'); % run Breach(Sys) to explore the successive counter-examples 
