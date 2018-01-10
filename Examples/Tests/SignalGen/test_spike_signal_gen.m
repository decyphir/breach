InitBreach
%%
%test_spike_signal_gen
input_gen = spike_signal_gen({'In1', 'In2'}, {'spline', 'linear'});
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.PrintParams();

%% Nominal trace
InputGen.Sim()
InputGen.PlotSignals()

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_base','In2_spike_base'}, [-1, 3]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_time','In2_spike_time'}, [5, 7]);
InputGen.Sim(); InputGen.PlotSignals();
%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_width','In2_spike_width'}, [12, 8]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_time','In2_spike_time', 'In1_spike_width','In2_spike_width'}, [5, 7, 12, 8]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_width','In2_spike_width'}, [0.1, 0.5]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_spike_amp','In2_spike_amp'}, [5, -5]);
InputGen.Sim();InputGen.PlotSignals();




