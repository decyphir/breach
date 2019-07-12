%% Breach Sets
%

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
clear; close all;
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

%% Specifying a domain (or type) for a parameter
% 
% We create an interface to the Abstract Fuel Control model as usual with a
% pulse input signal for pedal angle.
BrAFC.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'},...
                        [10 15; 0 61.1]);

%%
% We use SetDomain to specify that Pedal_Angle_pulse_period has to be
% integer. 
BrAFC.SetDomain('Pedal_Angle_pulse_period', 'int', [10 15]);

%% Sampling parameter sets with mixed domains
% The consequence of setting domains is that any subsequent sampling 
% restricts the parameters to belong to their respective domains. 
AFC_Grid= BrAFC.copy(); 

%Next we creates a 4x5 grid sampling and plot it.
AFC_Grid.GridSample([4 5]);  
figure; AFC_Grid.PlotParams({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'});  

%% More on domains and usage of function SampleDomain    

%% Domains  
% Breach supports integer and enum types  through so-called domains. 
p_int = 'Pedal_Angle_pulse_period';
p_enum = 'Pedal_Angle_pulse_amp';
p_double = 'Pedal_Angle_pulse_width';
p_double_bis = 'Pedal_Angle_base_value';
p_int_50 = 'Engine_Speed_u0';

%% 
% We set the corresponding domains:
BrAFC.SetDomain(p_int, 'int', [10 20]);
BrAFC.SetDomain(p_enum, 'enum', [10 50 80]);
BrAFC.SetDomain(p_double,'double', [0.1 1]);
BrAFC.SetDomain(p_double_bis,'double', [0 10]);
BrAFC.SetDomain(p_int_50, 'enum', 900:50:1100);  

%%
% There is also Boolean 'bool' type, which is equivalent to 'enum' with 
% values 0 or  1. 

%% SampleDomain examples
% SampleDomain is a function to sample domains. It is more generic,
% flexible and comprehensive than the legacy QuasiRandomSample and
% GridSample.

%%
% The general syntax is 

%% 
% B.SampleDomain(some_parameter, number_of_samples, method_of_sampling, behavior_combining_with_old_samples) 

%%
% Default method or sampling is random (uniform) and behavior for
% combining is 'replace'. E.g., 
BrAFC.SampleDomain(p_int,5);
BrAFC.GetParam(p_int) % get 5 random numbers between 10 and 20 

%% 
% Do it again:
BrAFC.SampleDomain(p_int,3); % get 3 random numbers between 10 and 20 
BrAFC.GetParam(p_int) % the 5 previous numbers were replaced

%% More examples
% Sample 4 values on a regular grid for p_double, and append the result with the 3 previous
% samples. Plotting makes it more clear. 
BrAFC.SampleDomain(p_double,4, 'grid', 'append');
figure; BrAFC.PlotParams({p_int, p_double})


%% More examples (ct'd)
% Combine previous samples with all values p_int_50
BrAFC.SampleDomain(p_int_50,'all',[],'combine');
figure; BrAFC.PlotParams({p_int, p_double, p_int_50});
view([45 45]);

%% More examples (ct'd)
% Sample all values of p_enum and combine with 4 random values of p_double
BrAFC.SampleDomain({p_enum, p_double},{'all', 10}, 'rand');
figure;BrAFC.PlotParams({p_enum, p_double}) 

%% More examples (ct'd)
% Simple 2d grid, 5x5
BrAFC.SampleDomain({p_double, p_double_bis} , 5, 'grid')
figure; BrAFC.PlotParams({p_double, p_double_bis})

%% More examples (ct'd)
% Simple 2d grid, 5x3
BrAFC.SampleDomain({p_double, p_double_bis} , [5 3], 'grid')
figure; BrAFC.PlotParams({p_double, p_double_bis})

%% More examples (ct'd)
% 100 random values in 2d domain. 
BrAFC.SampleDomain({p_double, p_double_bis} , 100)
figure; BrAFC.PlotParams({p_double, p_double_bis})





