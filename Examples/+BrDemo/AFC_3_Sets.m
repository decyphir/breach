%% Breach Sets
%

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
BrDemo.InitAFC;
BrAFC

%% BreachSets
% Each BreachSystem is also a BreachSet object, i.e., a collection of
% parameter values and traces. Some basic operations are available to
% generate new values, such as grid sampling and quasi-random sampling
% which we demonstrate next. 

% We first need to set the ranges of parameters that we want to vary.
AFC_Grid= BrAFC.copy(); 
AFC_Grid.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'},...
                        [10 15; 0 61.1]);
              
% Next we creates a 5x5 grid sampling and plot it.
AFC_Grid.GridSample([5 5]);  
AFC_Grid.PlotParams();  

%% Random Sampling of Parameters

% Again, we first set the ranges of parameters that we want to sample.
AFC_Rand= BrAFC.copy(); 
AFC_Rand.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp', 'Engine_Speed_u0'},...
                        [10 15; 0 61.1;900 1100]);
% Next we creates a random sampling with 50 samples and plot it.
AFC_Rand.QuasiRandomSample(50);  
figure; AFC_Rand.PlotParams();  
set(gca,'View', [45 45]);  

%% Simulation from Multiple Parameters (Grid)
AFC_Grid.Sim(0:.05:30);         
figure; AFC_Grid.PlotSignals({'Pedal_Angle','Engine_Speed', 'AF'}); % we plot only the input signals and AF

%% Simulation from Multiple Parameters (Random)
AFC_Rand.Sim(0:.05:30);         
figure; AFC_Rand.PlotSignals({'Pedal_Angle','Engine_Speed','AF'}); % we plot only the input signals and AF

%% Enumerating Corners

% Again, we first set the ranges of parameters that we want to sample.
AFC_Corners= BrAFC.copy(); 
AFC_Corners.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp', 'Engine_Speed_u0'},...
                        [10 15; 0 61.1;900 1100]);
% Next we obtain the corner parameters
AFC_Corners.CornerSample();  
figure;AFC_Corners.PlotParams();  
set(gca,'View', [45 45]);  

%% Simulations from Corners

AFC_Corners.Sim();
figure; AFC_Corners.PlotSignals({'Pedal_Angle','Engine_Speed','AF'});

