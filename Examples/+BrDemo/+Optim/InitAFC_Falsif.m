%% Initialize a number of falsification problems for fast benchmarking and testing

%% Model initialization 
BrDemo.InitAFC();

%% Setup an input space 
InputGen = fixed_cp_signal_gen({'Pedal_Angle','Engine_Speed'}, [3 1]);
BrAFC.SetInputGen(InputGen);
B = BrAFC.copy();
B.SetParamRanges({'Pedal_Angle_u1', 'Pedal_Angle_u2'}, [ 0 50 ]);
B.SetParamRanges({'Engine_Speed_u0'}, [ 900 1100 ] );
B.SetParam({'Pedal_Angle_u1', 'Pedal_Angle_u2', 'Engine_Speed_u0'}, [25; 25; 1000 ]);

STL_ReadFile('AFC_simple_spec.stl');

%% Easy 
phi005= set_params(AF_alw_ok, 'tol', 0.005);
AFC_falsif005 = FalsificationProblem(B, phi005);

%% Medium
phi01= set_params(AF_alw_ok, 'tol', 0.01);
AFC_falsif01 = FalsificationProblem(B, phi01);
phi02= set_params(AF_alw_ok, 'tol', 0.02);
AFC_falsif02 = FalsificationProblem(B, phi02);
phi03 = set_params(AF_alw_ok, 'tol', 0.03);
AFC_falsif03 = FalsificationProblem(B, phi03);

%% Hard to impossible
phi05= set_params(AF_alw_ok, 'tol', 0.05);
AFC_falsif05 = FalsificationProblem(B, phi05);

phi1= set_params(AF_alw_ok, 'tol', 0.1);
AFC_falsif1 = FalsificationProblem(B, phi1);

%% Requirement mining
Bsynth = BrAFC.copy();
Bsynth.AddSpec(AF_alw_ok);
Bsynth.SetParamRanges('tol', [0 0.1]);

synth_pb = ParamSynthProblem(Bsynth, AF_alw_ok, {'tol'},  [0 0.1]);
rq_mining_pb = ReqMiningProblem(AFC_falsif005, synth_pb);
