%% Create generator
input_gen = pulse_signal_gen({'In1', 'In2','In3'});
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.PrintParams();

%% Nominal trace
InputGen.Sim()
InputGen.PlotSignals()

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_base_value','In2_base_value','In3_base_value'}, [-1, 0, 3]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_pulse_period','In2_pulse_period','In3_pulse_period'}, [1, 2, 3]);
InputGen.Sim(); InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_pulse_width','In2_pulse_width','In3_pulse_width'}, [0.1, 0.5, 0.9]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_pulse_amp','In2_pulse_amp','In3_pulse_amp'}, [1, 5, 10]);
InputGen.Sim();InputGen.PlotSignals();

%% Varying parameters 
InputGen = BreachSignalGen({input_gen});  % generator with 3 pulse signals
InputGen.SetParam({'In1_pulse_delay','In2_pulse_delay','In3_pulse_delay'}, [-0.2, 0.3, 2]);
InputGen.Sim();InputGen.PlotSignals();


