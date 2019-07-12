%% Signal Temporal Logic (STL) Specifications 

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
clear; close all;
BrDemo.InitAFC;
BrAFC

%% Writing a Simple STL Specification 
% First we define a predicate stating that AF is above 1% of AFref
AF_not_ok = STL_Formula('AF_ok', 'abs(AF[t]- AFref[t]) < 0.01*14.7')

%%
% The first argument 'AF_ok' is an id for the formula needed so that
% it can be referenced as a sub-formula. 

%%
% The second argument is a string describing the formula. The syntax AF[t] refers to the value of  
% signal AF at a time t. If we evaluate this formula, t will be
% instantiated with 0. We need temporal operators to examine other time
% instants. 

%%
% For example, we define next a formula stating "some time in the future, AF_ok is true". 
AF_ev_ok = STL_Formula('AF_ev_ok', 'ev (AF_ok)')

%%
% ev is a shorthand for 'eventually'. Other temporal operators are 'alw' or 'always' and 'until'.  

%% Checking a Simple Specification on a Simulation
% Temporal operators can specify a time range to operate on. For our
% system, AF does not need to be checked before 10s, and we simulate until
% 30s so we check the following formula:
AF_alw_ok = STL_Formula('AF_alw_ok', 'alw_[10,30] (AF_ok)') % alw shorthand for in 'always'

%%
% Then we check our formula on a simulation with nominal parameters:
AFC_w_Specs= BrAFC.copy();
Time = 0:.05:30;
AFC_w_Specs.Sim(Time);
AFC_w_Specs.CheckSpec(AF_alw_ok)

%% 
% A positive results means that the formula is satisfied, i.e., the system does not 
% overshoot. 

%% Checking a Simple Specification on a Simulation (Plot)
% We can plot the satisfaction function with 
figure; AFC_w_Specs.PlotRobustSat(AF_alw_ok);

%% Checking Another Formula
% The reason we are not interested in AF between 0 and 10s is because the 
% controller is not in a mode where it tries to regulate it at this time.
% We can implement this explicitly using the controller_mode signal:
AF_alw_ok2 = STL_Formula('AF_alw_ok2', 'alw (controller_mode[t]==0 => AF_ok)')

%%
% Then check the new formula on the simulation we performed already:
AFC_w_Specs.CheckSpec(AF_alw_ok2)

%% 
% This formula is not satisfied - something might have gone wrong
% somewhere. 

%% Plotting Satisfaction for Debugging Formula 
% At time t=0, controller_mode is briefly 0, which causes the negative results, but this is only an initialization glitch.  
figure; AFC_w_Specs.PlotRobustSat(AF_alw_ok2);

%% Reading a formula from an STL File
% STL formulas can be defined in a file. This makes it much easier to write
% complex properties. Also easier to defined parametric STL formulas, e.g.: 
%
type AFC_simple_spec.stl

%%  
% STL files are loaded using the following command:

STL_ReadFile('AFC_simple_spec.stl'); % this loads all formulas and sub-formulas defined in the file
AFC_w_Specs.CheckSpec(AF_alw_ok)


%% Plotting Satisfaction    
% We can plot the satisfaction function again to interpret the result.  
figure; AFC_w_Specs.PlotRobustSat(AF_alw_ok);

%% Specifying Depth for Satisfaction Plots
% PlotRobustSat decomposes the formula into subformulas. We can specify the
% depth of decomposition, e.g.,

figure; AFC_w_Specs.PlotRobustSat(AF_alw_ok,3); % shows satisfaction of top formula and 2 subformulas.

%% Formula Parameters with set_params/get_params
% Parameters for formulas can be accessed and changed using get_params and set_params functions: 
get_params(AF_alw_ok)

%%
%
AF_alw_ok2 = set_params(AF_alw_ok, {'tol'}, [1e-3]) % copy formula with tolerance changed from 1% to 0.1%
get_params(AF_alw_ok2)

%% Formula Parameters with GetParam/SetParam
% Another way to access and modify formula parameters is by defining them
% using the method SetParamSpec of class BreachSystem. This
% actually overrides the values defined for the STL_Formula object.

AFC_w_Specs.SetParamSpec('tol', 1e-3);
AFC_w_Specs.CheckSpec(AF_alw_ok)

%%
% The advantage of the second method is that we can explore (sample) the formula
% parameter space the same way we did with system parameters, as
% demonstrated next.   

%% Checking Multiple Formula Parameters

AFC_w_mult_Specs= AFC_w_Specs.copy();
% we define a range for parameter tol and and grid sample it with 10 values 
AFC_w_mult_Specs.SetParamRanges('tol', [0 1e-2]); 
AFC_w_mult_Specs.GridSample(10); AFC_w_mult_Specs.GetParam('tol')

%%
%
AFC_w_mult_Specs.CheckSpec(AF_alw_ok)

%% Checking a Specification for Multiple System Parameters

% Create new interface object with some ranges for pedal angle inputs and sample randomly
AFC_Rand_w_Specs = BrAFC.copy();
AFC_Rand_w_Specs.SetParamRanges({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [0 20; 10 15; 0 40]);
AFC_Rand_w_Specs.QuasiRandomSample(10);
AFC_Rand_w_Specs.Sim(Time);

%% Set specification parameter
AFC_Rand_w_Specs.SetParam({'t_start', 't_end', 'tol'}, [10 30 .003], 'spec'); 
% Checks Specification for all simulations
AFC_Rand_w_Specs.CheckSpec(AF_alw_ok)
%%
% Note that for some simulations, the specification is falsified, meaning the system
% overshoots. 

%% Checking a Specification for Multiple System Parameters
% To find out which parameter value lead to a positive or negative
% satisfaction, we can do the following:

[Brpos, Brneg] = AFC_Rand_w_Specs.FilterSpec(AF_alw_ok);
good_param = Brpos.GetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'})
bad_param = Brneg.GetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'})


%% Monitoring a Formula on a Trace (1)
% To monitor STL formulas on existing traces, one can use the
% BreachTraceSystem class. 

% create a trace with signals x and y
time = 0:.05:10; x = cos(time); y = sin(time);
trace = [time' x' y']; % trace is in column format, first column is time
BrTrace = BreachTraceSystem({'x','y'}, trace); 
figure; BrTrace.PlotSignals();

%% Monitoring a Formula on a Trace (2)
% Checks (plots) some formula on imported trace:
figure; BrTrace.PlotRobustSat('alw (x[t] > 0) or alw (y[t]>0)');





