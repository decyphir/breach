%% Breach demo: analysis of an Automatic transmission controller

%% Interface Automatic Transmission model with Breach 
B = BreachSimulinkSystem('Autotrans_shift');
B.PrintAll % print available signals and parameters 

%% Running one simulation
B.SetTime(0:.01:30); B.SetParam({'throttle_u0'}, 100);
B.Sim(); B.PlotSignals({'throttle', 'RPM', 'speed', 'gear'});

%% Checks a property: The speed is never below 30 while in gear3 
STL_ReadFile('Autotrans_spec.stl');
B.PlotRobustSat(gear3_and_speed_low)

%% Describes and generate driving scenarios 
% We create an input generator that will alternates between acceleration and braking 
sg = var_step_signal_gen({'throttle', 'brake'}, 5);
B.SetInputGen(sg);
            
% We assign ranges for duration and amplitude of each input:
B.SetParamRanges({'dt_u0', 'dt_u1', 'dt_u2', 'dt_u3'}, ...
                  [.1 10  ;  .1 10;    0.1 10;    0.1 10]);
B.SetParamRanges({'throttle_u0','brake_u1', 'throttle_u2', 'brake_u3'}, ... 
                  [0 100;        0 325;      0 100;         0 325]);

% We don't specify a range for brake_u0 so that it remains constant equal
% to 0 (by default). Same for throttle_u1, etc.
B.QuasiRandomSample(10); B.Sim();

%% Plot multiple simulations result
B.PlotSignals({'throttle', 'brake','RPM', 'speed', 'gear'});

%% Check property visually 
figure;
B.PlotSigPortrait({'gear', 'speed'})   

%% Checks property by monitoring 
B.CheckSpec(never_gear3_and_speed_low);
B.PrintSpecs

%% Falsify property
B.ResetSampling(); % remove the 10 samples and traces, keep parameter ranges
falsif_pb = FalsificationProblem(B, never_gear3_and_speed_low);
falsif_pb.max_obj_eval = 1000;
falsif_pb.max_time = 180; % give the solver three minutes to falsify the property
falsif_pb.solve();

%% Examine counter-example
BFalse = falsif_pb.BrSet_False;
BFalse.PlotSigPortrait({'gear','speed'});

%% Fix Specification
param_pb = ParamSynthProblem(BFalse, never_gear3_and_speed_low, 'v_low', [0 30]);
param_pb.solve();

%% Requirement mining: Iterate 
mining_pb = ReqMiningProblem(param_pb, falsif_pb);
mining_pb.solve();

%% Visualizing final coverage 
mining_pb.falsif_pb.BrSet_Logged.PlotSigPortrait({'gear','speed'})

