
%% 
Br = BreachSimulinkSystem('Autotrans_shift');
Br.SetTime(0:.01:50);

%% load formula
formulas = STL_ReadFile('spec.stl');
phi_template = phi1;
      
%% Property parameters 
prop_params.names = {'vmax', 'rpm_max'};
prop_params.ranges = [0 200   ;...  % for vmax
                      0 6000 ];     % for rmp_max
  
%% Input parameters
input_params.names = {'throttle_u0' ... ,
                     'brake_u0'... 
                    };
input_params.ranges = [0 100; ...
                     0 325; ...   
];

%% Compute some traces 
Br.SetParamRanges(input_params.names, input_params.ranges);
Br.QuasiRandomSample(10);
Br.Sim(); 

%% Param synthesis problem 

synth_problem = ParamSynthProblem(Br, phi1, prop_params.names, prop_params.ranges);
synth_problem.setup_solver('fmincon')
synth_problem.solve();

%% Get param synthesis results

