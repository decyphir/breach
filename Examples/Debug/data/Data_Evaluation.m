% Set things up
clear;
InitBreach;

% Load the data
M = csvread('Data_1.csv');
x = M(:,[1:2]);
g = M(:,[1,3]);
t = transpose(M(:,1));

% Interface workspace data with a Breach object
workspace_data = from_workspace_signal_gen({'x','g'});
B = BreachSignalGen({workspace_data});
B.Sim(t);
 
% Load STL file and associate with the Breach object
STL_ReadFile('Data_Spec.stl');

% Compute qualitative and quantitative semantics
close all;
warning('off','STL_Eval:Inf_or_Nan'); % We use Inf a lot and some NaNs also 
                                      % seem to be generated in this case

% Robustness as computed by default in Breach
figure('Name', 'Standard Robsutness');
B.PlotRobustSat(phi);

%% The two most interesting IO-robustness notions:
% Relative output robustness: for a fixed input, how much to change the
% output so as to change the satisfaction status of the formula
% This can serve as a system robsutness measure
figure('Name', 'Relative Output Robustness');
B.PlotIORobustSat(phi, 'out', 'rel');

% Absolute input robustness: how much to change the
% input so as to change the satisfaction status of the formula for all
% possible outputs
% This can serve as a formula vacuity measure
figure('Name', 'Absolute Input Robustness');
B.PlotIORobustSat(phi, 'in', 'abs');

%% Two other notions of robustness that may or may not be interesting:
% Relative input robustness: for a fixed output, how much to change the
% input so as to change the satisfaction status of the formula
figure('Name', 'Relative Input Robustness');
B.PlotIORobustSat(phi, 'in', 'rel');

% Absolute output robustness: how much to change the
% output so as to change the satisfaction status of the formula for all
% possible inputs
figure('Name', 'Absolute Output Robustness');
B.PlotIORobustSat(phi, 'out', 'abs');