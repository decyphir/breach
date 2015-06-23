
%% Name of the Simulink model 
mdl = 'helicopter';

%% Parameters needed to run the model 
K = 10;
i = 0; 
a = 1;

%% Interface with Breach 
Sys = CreateSimulinkSystem(mdl)

%% Default 
Sys.tspan = 0:.01:4;

%% Default value for input signal
Sys = SetParam(Sys, 'psi_u0', 10);