InitBreach;

%% load formula
formulas = QMITL_ReadFile('spec.stl');
phi_template = phi1;

%% Create system and input strategy 
mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
Sys.tspan = 0:.01:50;
      
%% Property parameters 
prop_opt.params = {'vmax', 'rpm_max'};
prop_opt.monotony   = [1 1];
prop_opt.p_interval = [0 200   ;...  % for vmax
                       0 6000 ];     % for rmp_max
prop_opt.p_tol      = [1 1];
  
%% System and falsification parameters
falsif_opt.params = {'throttle_u0' ... ,
                     'brake_u0'... 
                    };
falsif_opt.ranges = [0 100; ...
                     0 325; ...   
];
falsif_opt.nb_init = 10;
falsif_opt.nb_iter = 1000;
falsif_opt.nb_max_call = 1000;

%% Max number of mining iterations
iter_max= 100;
%[p, rob, Pr] = ReqMining(Sys, phi_template, falsif_opt, prop_opt, iter_max);