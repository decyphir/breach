%% Advanced Analysis

%% Initialization
BrDemo.InitAFC;

%%
% 

% Load some properties.
STL_ReadFile('AFC_simple_spec.stl');

%% Plotting Satisfaction Maps
% We can visualize the satisfaction of a formula when varying
% parameters. For example:
BrMap = BrAFC.copy(); 
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp'}, [10 80]);

%% Plotting Satisfaction Maps (ct'd)
% The Robust Map can be plotted as a surface plot, varying two parameters:
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp', 'Pedal_Angle_pulse_period'}, [0 80; 10 20]);

%% Plotting Satisfaction Maps (ct'd)
% We can zoom in to get more details. The last argument of PlotRobustMap
% can be used to get a more refined grid too (here 15x15).
BrMap.PlotRobustMap(AF_alw_ok, {'Pedal_Angle_pulse_amp', 'Pedal_Angle_pulse_period'}, [40 60; 15 20], [15 15]); 

%% Sensitivity Analysis
BrSensi = BrAFC.copy();
params = {'AF_sensor_tol','MAF_sensor_tol','fuel_inj_tol','kappa_tol','tau_ww_tol','pump_tol'};
ranges = [ 0.99 1.01; 0.99 1.01; 0.99 1.01;0.99 1.01; 0.99 1.01;0.99 1.01];
BrSensi.SensiSpec(AF_alw_ok, params, ranges);
