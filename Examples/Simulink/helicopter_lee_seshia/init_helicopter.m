%% Creating an interface with Breach using CreateSimulinkSystem 

%%
% Name of the Simulink model:
mdl = 'helicopter';

%%
% Next, we define the parameters needed to run the model in the workspace.
% These parameters will be detected by Breach, and included in its 
% interface with the system. 
K = 10;
i = 0; 
a = 1;

%%
% By default, Breach will look for inputs, outputs, logged signals and 
% scoped signals (i.e. signals that are inputs to a scope block in the 
% model. Breach interface will then allow to change parameters and constant 
% input signals, simulate the model and collect and observe the resulting 
% trace. 
Sys = CreateSimulinkSystem(mdl);

%%
% We set the default simulation time:
Sys.tspan = 0:.01:4;

%%
% We set the default value for the (constant) input signal: 
Sys = SetParam(Sys, 'psi_u0', 10);
