BrDemo.InitAFC
STL_ReadFile('AFC_simple_spec.stl');

BrMap = BrAFC.copy(); 

%%
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp'}, [10 80]);

%%
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp', 'Pedal_Angle_pulse_period'}, [0 80; 10 20]);

%%
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp', 'Pedal_Angle_pulse_period'}, [40 60; 15 20]);
