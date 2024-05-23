BP2param
model_name = 'BP2_in';
IO_signal_names = {'In1','OutSat'};
warning('OFF', 'Simulink:blocks:NonPositiveIntegerValueSeed');
BP2 = BreachSimulinkSystem(model_name, {}, [],  IO_signal_names);
