% This is an example script of how to use Breach with the Turbo solver. 

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
BrDemo.InitAFC();
BrAFC

%% Creating a Falsification Problem
% Given a requirement R and some parameter range, we want to find a
% parameter value for which the system violates R.

%% 
% First we create the parameter search domain and load specifications. 
AFC_Falsify = BrAFC.copy();
AFC_Falsify.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 12; 10 40]);
STL_ReadFile('AFC_simple_spec.stl');
req = BreachRequirement(AF_alw_ok);

savingExampleName = fopen('currentExample.txt','wt');
fprintf(savingExampleName, 'demo_example');
fclose(savingExampleName);
%% 
% Then we create the falsification problem and solve it.
falsif_pb = FalsificationProblem(AFC_Falsify, req);
falsif_pb.max_obj_eval = 23;

%% Calling turbo
xHist = [10 12 10
        10 40 40];
fHist = [0.01992 0.11524 0.014005];
falsif_pb.setup_turbo('start_sample', xHist, 'start_function_values', fHist);
falsif_pb.solve();
