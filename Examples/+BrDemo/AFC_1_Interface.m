%% Demo 1: Interfacing a Simulink model
% Wherein we demonstrate how to interface a Simulink model with Breach. 

%% Interfacing a Simulink model with Breach  
% Breach can interact with a Simulink model by 
%
% * Changing parameters defined in the base workspace
% * Generating input signals  
% * Running simulations and collecting signals for analysis
%  
% We use the model AbstractFuelControl as running example.

%%
% First, we initialize Breach and instantiate parameters for the model to run properly: 
clear;close all;
InitBreach;
fuel_inj_tol = 1.0; 
MAF_sensor_tol = 1.0;
AF_sensor_tol = 1.0; 
pump_tol = 1.;
kappa_tol=1; 
tau_ww_tol=1;
fault_time=50;
kp = 0.04;
ki = 0.14;

%%
% They will be discovered automatically and included in the interface. Good
% practice is to create a separate script initializing parameters of
% concern, see for example BrDemo.InitAFCparams.

%% Creating an Interface Object (1)
%
% Running the following will create a BreachSystem interface with the model:
mdl = 'AbstractFuelControl';
warning('off', 'Simulink:LoadSave:EncodingMismatch') %avoid warning encoding conflict windows vs iso
BrAFC = BreachSimulinkSystem(mdl)

%% 
% Note that Breach created a copy of the model in AbstractFuelControl_breach

%%
% We can check which parameters are interfaced via the following command: 
BrAFC.PrintParams();


%% Creating an Interface Object (2)
% 
% Sometimes a model contains many tunable parameters. In that case, we can
% explicitly specify those we want to tune, e.g., some PI parameters:

BrAFC_less_params = BreachSimulinkSystem(mdl, {'ki', 'kp'}, [0.14 0.04]);
BrAFC_less_params.PrintParams();


%% 
% Similarly,we may want to have only a few signals among the ones logged
% to be visible:

BrAFC_less_signals = BreachSimulinkSystem(mdl, 'all', [], ... % 'all' means detects all visible parameters
                                             {'Pedal_Angle','Engine_Speed', 'AF'}); 
BrAFC_less_signals.PrintAll;


%% Interfacing Signals
% By default, Breach collects the following signals:
%
% * Signals attached to input ports 
% * Signals attached to output ports
% * Logged signals 
% 

%% 
% We can check which signals are interfaced via the following command: 
BrAFC.PrintSignals();

%%
% Note that signal and parameter names in the model should remain simple. 
% Preferably valid variable id names, i.e., using only underscores and 
% alphanumerical characters (a-z, A-Z, 0-9 or '_').  


%% Input Signals
% Input signals of different types can be generated. They are defined by 
% parameters of the form input1_u0, input1_u1 etc.

%% 
% By default, input signals are simply constant and we have two inputs,  
% which are Engine_Speed and Pedal\_Angle. Next, 
% we set their constant values: 
BrAFC.SetParam('Engine_Speed_u0',1000)
BrAFC.SetParam('Pedal_Angle_u0',30)
BrAFC.PrintParams();

%%
% We can now run simulations.  

%% Simulation with Nominal Parameters
%

AFC_Nominal= BrAFC.copy(); % Creates a copy of the interface object 
AFC_Nominal.Sim(0:.05:30); % Run Simulink simulation until time=30s 
AFC_Nominal.PlotSignals(); % plots the signals collected after the simulation


%% Plotting Nominal Simulation
% The plotting methods accepts a number of arguments. The first one and
% most useful simply allows to select which signals to plot. 

figure; AFC_Nominal.PlotSignals({'Pedal_Angle','Engine_Speed','cyl_air', 'cyl_fuel', 'AF'}); 


