%% Create generator
input_gen = multi_pulse_signal_gen({'In1', 'In2','In3'});
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.PrintParams();

%% Nominal trace
figure;
InputGen.Sim()
InputGen.PlotSignals()

%% Varying parameters 
figure;
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'period','delay','In2_pulse_rel_delay','In3_pulse_rel_delay'}, [1.5, .5, 0.25, 0.5]);
InputGen.Sim();InputGen.PlotSignals();
