
InitBreach;
%% load formula
formulas = STL_ReadFile('spec.stl');
phi_template = phi_stay2;

%% Create system and input strategy 
mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'VarStep2');
Sys.tspan = 0:.01:50;

% Accelerate for 20 s at 100%
Sys = SetParam(Sys, 'throttle_dt0', 20);
Sys = SetParam(Sys, 'throttle_u0', 100);
Sys = SetParam(Sys, 'brake_dt0', 20);
Sys = SetParam(Sys, 'brake_u0', 0);
     
% Brake for 20s after
Sys = SetParam(Sys, 'throttle_dt1', 20);
Sys = SetParam(Sys, 'throttle_u1',   0);
Sys = SetParam(Sys, 'brake_dt0', 20);
Sys = SetParam(Sys, 'brake_u0', 325);

%% Property parameters 
prop_opt = []; 
prop_opt.params = { 't1'};
prop_opt.monotony = [-1];
prop_opt.p_tol = [.02];
prop_opt.ranges = [0 40]
               
%% System and falsification parameters
falsif_opt = [];
falsif_opt.params = {'throttle_dt0', 'throttle_u0','brake_dt0','brake_u1'};
falsif_opt.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 20   ;
                      0 325 ;
                       ];
 
falsif_opt.nb_init =  1000;
falsif_opt.nb_iter = 1000;  
falsif_opt.nb_max_call = 10000;

iter_max= 100;
[p,rob,Pr_stay_gear2] = ReqMining(Sys, phi_template, falsif_opt, prop_opt,iter_max)
Psave(Sys, 'Pr_stay_gear2'); % run Breach(Sys) to explore the successive counter-examples 
