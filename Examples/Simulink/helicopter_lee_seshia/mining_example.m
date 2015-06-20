clear 
InitBreach;

%% Init system
init_helicopter;

%% Set default input 
Sys = SetParam(Sys, 'psi_u0', 10);

%% Compute nominal behavior
P = CreateParamSet(Sys);
P = ComputeTraj(Sys,P);
SplotVar(P);
pause
close

%% load properties
QMITL_ReadFile('specs.stl');

%% Parameter synthesis for nominal trajectory

% Property parameters 
prop_opt.params = {'tau'};
prop_opt.monotony   = [1];
prop_opt.ranges = [0 10];
prop_opt.p_tol = [.01];

[p, rob, P] = GetPropParamBin(Sys, phi, P, prop_opt);

figure;
SplotSat(Sys, P, phi, 3);
%title('Candidate specification for nominal trace');
pause
close

%% System and falsification for nominal parameter
falsif_opt.params = {'K' ... ,
                    };

falsif_opt.ranges = [9 11; ...
];

falsif_opt.nb_init = 10;
falsif_opt.nb_iter = 10;
falsif_opt.nb_max_call = 100;

prop_opt.values= p;
Pfalse = Falsify(Sys, phi, falsif_opt, prop_opt);

figure;
SplotSat(Sys, Pfalse, phi, 3);
%title('Falsifying trace');
pause
close 

%% Max number of mining iterations
iter_max= 10;
[p, rob, Pr] = ReqMining(Sys, phi, falsif_opt, prop_opt, iter_max);
Psave(Sys, 'Pr');

%% Ground truth
Pground = SetParam(P,'K', 9); 
Pground = ComputeTraj(Sys,Pground);

[p_ground, rob, Pground] = GetPropParamBin(Sys, phi, Pground, prop_opt);






