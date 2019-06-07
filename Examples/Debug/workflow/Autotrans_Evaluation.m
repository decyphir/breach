%% Interface Automatic Transmission model with Breach 
B = BreachSimulinkSystem('Autotrans_Model');
B.PrintAll % print available signals and parameters 

%% Running one simulation
B.SetTime(0:.01:30); B.SetParam({'throttle_u0'}, 100);
B.Sim(); B.PlotSignals({'throttle', 'RPM', 'speed', 'gear'});

%% Checks a property: The speed is never below 30 while in gear3 
STL_ReadFile('Autotrans_Spec.stl');
warning('off','STL_Eval:Inf_or_Nan');
close all;
figure; B.PlotIORobustSat(gear3_and_speed_low, 'in', 'rel');
figure; B.PlotIORobustSat(gear3_and_speed_low, 'in', 'abs');
figure; B.PlotIORobustSat(gear3_and_speed_low, 'out', 'rel');
figure; B.PlotIORobustSat(gear3_and_speed_low, 'out', 'abs');