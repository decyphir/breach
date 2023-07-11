
addpath('/Users/thaodang/Metaheuristics/src')
addpath('/Users/thaodang/Metaheuristics/wordgen')
addpath('/Users/thaodang/Metaheuristics/breach-dev')
addpath('.')    
InitBreach('/Users/thaodang/Metaheuristics/breach-dev',true); % forces initialization from folder in Metaheuristics


%% TA based signal generator setup 
%

num_evt=5

init_TA_signal_gen;


%%  Checking reachable labels
STL_ReadFile('Autotrans_req.stl');

%% 

S2 = S0.copy();
mdl = 'Autotrans_wordgen';
%mdl = 'autotrans_mod04';
Ba = BreachSimulinkSystem(mdl);
Ba.SetTime(time);
Ba.SetInputGen(S2);
Ba.SetParamRanges([pevts pbranching], [0 1]);
Ba.Sim();


%% 
R = BreachRequirement(never_gear3_and_speed_low);
falsif_pb = FalsificationProblem(Ba, R);
falsif_pb.solver_options.num_corners = 100;
falsif_pb.solver_options.num_quasi_rand_samples = 100;
falsif_pb.max_obj_eval = 1000;
%falsif_pb.SetupDiskCaching();

%%
falsif_pb.solve();

%% Initial Counter-example
BFalse = falsif_pb.BrSet_False;
%BFalse = falsif_pb.GetFalse();;

%% Fix Specification
param_pb = ParamSynthProblem(BFalse, never_gear3_and_speed_low, 'v_low', [0 30]);
param_pb.solver_options.monotony = -1;
param_pb.solve();


%% Requirement mining: Iterate 
%mining_pb = ReqMiningProblem(param_pb, falsif_pb);
%mining_pb.solve();
