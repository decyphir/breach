%% Online Monitoring 
% 

%% Initialization
% This demo uses a modified version of the AbstractFuelControl model to
% illustrate the use of online monitoring blocks from the Breach Library.
BrDemo.InitAFC_Online()

%%
% Breach Library is accessible via Simulink Library Browser, or by opening 
% slstlib.slx. AFC_Online contains one 'STL Monitor' block and one 'STL
% stops when false' block. We configure the first one as follows.

max_rob = 0.5;            % estimate of maximum robust satisfaction of a formula (default to +inf)
sig_names = 'AF,AFref';   % declares signal names used in the specification

% defines a simple overshoot property - note that STL formula are given as
% string, STL_Formula object are not supported. 
phi_st = 'alw_[10, 30] ((abs(AF[t]-AFref[t]) > 0.05) => (ev_[0, 1] (abs(AF[t]-AFref[t]) < 0.05)))';

%%
% Set parameters for one simulation and run. 
BrAFC_Online.SetParam({'max_rob','Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [max_rob, 12, 50]);
BrAFC_Online.Sim(0:.1:40) 

%% Plotting interval robustness 
% At every time satisfaction value is guaranteed to be within upper and lower
% robustness.
hold on; BrAFC_Online.PlotExpr({'abs(AF[t]-AFref[t])','rho_up[t]', 'rho_low[t]'}, 'LineWidth',2);
grid on;set(gca, 'YLim',[-0.6, 0.6]);
legend({'abs(AF-AFref))', 'Upper robustness','Lower robustness'});



%% Stopping simulation as soon as possible
% As soon as upper robustness becomes negative, simulation is stopped by
% the 'Stop when false block'. We define its specification as follows:
phi_stop = '(abs(AF[t]-AFref[t]) > 0.05) => (ev_[0, 1] (abs(AF[t]-AFref[t]) < 0.05))';
BrAFC_Online.ResetSimulations();
BrAFC_Online.Sim(0:.1:40) 
%%
% Note that the 'Stop when false block' checks its formula at every time,
% whereas the STL Monitor block checks it at time 0. 

%% Stopping simulation as soon as possible - plot
%
figure; hold on; BrAFC_Online.PlotExpr({'AF[t]-AFref[t]','rho_up[t]', 'rho_low[t]'}, 'LineWidth',2);
grid on;set(gca, 'YLim',[-0.6, 0.6]);
legend({'abs(AF-AFref))', 'Upper robustness','Lower robustness'});

%% Online Monitoring Limitations 
%

%% 
% The online monitoring blocks do not support full parametric STL as
% offline monitoring. Below are the main limitations:
%
% * Atomic predicates only support operations: +,-,*,/,abs 
% * No parameters in predicate or time intervals  
% * No retiming for signals. E.g., x[t+0.1] is not supported
% 
% Also note that specifications are provided as strings, not as STL_Formula
% objects. In the dialog box, do not forget quotes, e.g., 'x1[t]>0'. 


