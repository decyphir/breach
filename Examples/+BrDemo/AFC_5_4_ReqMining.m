%% Requirement Mining Problems

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
BrDemo.InitAFC;
BrAFC

%%
% Load AFC specs.
STL_ReadFile('AFC_simple_spec.stl');

%% Requirement Mining
% We falsified AF_alw_ok, then updated the tol parameter to make it true again.

%%
% A natural next step would be falsify the property with new
% tolerance, and iterate. We can automate this process using a ReqMiningProblem. 
AFC_ReqMining = BrAFC.copy();

% Input parameters names and ranges for falsification
input_params = struct('names',{{'Pedal_Angle_pulse_period','Pedal_Angle_pulse_amp'}}, ...                    
                      'ranges', [10 20; 10 60]);

% Property with parameters names and ranges for paramter synthesis
phi = AF_alw_ok;
prop_params = struct('names', {{'tol'}}, 'ranges', [0 0.1]); 

%% Requirement Mining (ct'd)
% A ReqMiningProblem is an object combining a FalsificationProblem object
% and a ParamSynthesis object.
mine_pb = ReqMiningProblem(AFC_ReqMining, phi, input_params, prop_params);

%%
% We can change options for both problems as usual. E.g., specify
% monotonicity:
mine_pb.synth_pb.solver_options.monotony = 1;

%% Requirement Mining - Solving
mine_pb.solve();

%% Requirement Mining
% The mining ends when the falsifier fails. 