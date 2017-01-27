%% Parameter Synthesis Problems

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
BrDemo.InitAFC;
BrAFC

%%
% 

% Load some properties.
STL_ReadFile('AFC_simple_spec.stl');


%% Parameter Synthesis
% We falsified the property that AF stays within tol=0.01 of
% AF_ref. We might also ask: for what minimum value of tol is this property true?

%%
% This can be solved by solving a parameter synthesis problem:

AFC_ParamSynth = BrAFC.copy(); AFC_ParamSynth.Sim();
synth_pb = ParamSynthProblem(AFC_ParamSynth, AF_alw_ok, {'tol'},  [0 0.1]);
synth_pb.solve();

%%
% The default solver guessed that the problem was monotonic in 'tol'.
% Monotonicity is actually required by this solver. 

%%
% If monotonicity cannot be guessed it can be set up manually:
synth_pb.solver_options.monotony = 1; % +1/-1 resp. increasing/decreasing

%%
% If the problem is not monotonic, another solver has to be chosen.

%% Parameter Synthesis - Plot
AFC_ParamSynth.PlotSignals({'AF'}, [], {'LineWidth', 2});
set(gca, 'XLim', [10 40], 'FontSize',14, 'LineWidth', 2);

% plotting the new tolerance region
tol_best  = synth_pb.x_best;
plot([0 41], (1+tol_best)*[14.7 14.7],'r'); plot([0 41], (1-tol_best)*[14.7 14.7],'r');

%% Parameter Synthesis For Multiple Traces
% Parameter synthesis can be solved for multiple traces. 

AFC_ParamSynth.SetParamRanges({'Pedal_Angle_pulse_period','Pedal_Angle_pulse_amp'},[10 20;10 60]);
AFC_ParamSynth.QuasiRandomSample(10);
AFC_ParamSynth.Sim();
synth_pb = ParamSynthProblem(AFC_ParamSynth, AF_alw_ok, {'tol'},  [0 0.1]);
synth_pb.solve();

%% Parameter Synthesis For Multiple Traces - Plot
AFC_ParamSynth.PlotSignals({'AF'});
set(gca, 'XLim', [10 40]);
% plotting the new tolerance region
tol_best  = synth_pb.x_best;
plot([0 41], (1+tol_best)*[14.7 14.7],'r'); plot([0 41], (1-tol_best)*[14.7 14.7],'r');



