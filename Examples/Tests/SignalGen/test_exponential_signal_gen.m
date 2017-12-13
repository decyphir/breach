InitBreach
%%
%test_spike_signal_gen
input_gen = exponential_signal_gen({'In1', 'In2'});
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.PrintParams();

%% Nominal trace
InputGen.Sim()
InputGen.PlotSignals()

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_exp_base','In2_exp_base'}, [-1, 3]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_exp_amp','In2_exp_amp'}, [5, 20]);
InputGen.Sim(); InputGen.PlotSignals();
%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_exp_rate','In2_exp_rate'}, [0.2, 0.01]);
InputGen.Sim();InputGen.PlotSignals();
%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_exp_rate','In2_exp_rate'}, [-0.5, -0.02]);
InputGen.Sim();InputGen.PlotSignals();





