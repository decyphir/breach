InitBreach
%%
%test_spike_signal_gen
input_gen = sinusoid_signal_gen({'In1', 'In2'});
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.PrintParams();

%% Nominal trace
InputGen.Sim()
InputGen.PlotSignals()

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_sin_base','In2_sin_base'}, [-1, 3]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_sin_amp','In2_sin_amp'}, [5, -7]);
InputGen.Sim(); InputGen.PlotSignals();
%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_sin_decay','In2_sin_decay'}, [0.1, -0.1]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_sin_freq','In2_sin_freq'}, [3, 1.5]);
InputGen.Sim();InputGen.PlotSignals();




