%% Creates interface and run one simulation
Br = BreachSimulinkSystem('Autotrans_shift');
Br.SetTime(0:.01:50);

% Simulate once
Br.Sim();

%% Defines the Req mining problem
% Input parameters
input_params.names = {'throttle_u0' ... ,
                     'brake_u0'... 
                    };
input_params.ranges = [0 100; ...
                     0 325; ...   
];

% Template property
STL_ReadFile('spec.stl');
phi_template = phi1;
      
% Property parameters 
prop_params.names = {'vmax', 'rpm_max'};
prop_params.ranges = [0 200   ;...  % for vmax
                      0 6000 ];     % for rmp_max

                             
mine_phi1 = ReqMiningProblem(Br, phi_template, input_params, prop_params);

%% Solve it
mine_phi1.solve();
