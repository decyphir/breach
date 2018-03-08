%% Initialize a number of falsification problems for fast benchmarking and testing

%% Model initialization 
BrDemo.InitAFC();

%% Setup an input space 
InputGen = fixed_cp_signal_gen({'Pedal_Angle','Engine_Speed'}, [3 1]);
BrAFC.SetInputGen(InputGen);
B = BrAFC.copy();
B.SetParamRanges({'Pedal_Angle_u1', 'Pedal_Angle_u2'}, [ 0 50 ]);
B.SetParamRanges({'Engine_Speed_u0'}, [ 900 1100 ] );
STL_ReadFile('AFC_simple_spec.stl');

%% Easy 
phi005= set_params(AF_alw_ok, 'tol', 0.005);
AFC_falsif005 = FalsificationProblem(B, phi005);

%% Medium
phi05= set_params(AF_alw_ok, 'tol', 0.05);
AFC_falsif05 = FalsificationProblem(B, phi05);

%% Hard to impossible
phi1= set_params(AF_alw_ok, 'tol', 0.1);
AFC_falsif1 = FalsificationProblem(B, phi1);

%% Requirement mining
Bsynth = BrAFC.copy();
Bsynth.AddSpec(AF_alw_ok)
Bsynth.SetParamRanges('tol', [0 0.1]);

synth_pb = ParamSynthProblem(Bsynth, AF_alw_ok, {'tol'},  [0 0.1]);
rq_mining_pb = ReqMiningProblem(AFC_falsif005, synth_pb);
