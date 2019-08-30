%% Falsification Problems

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
clear; close all;
BrDemo.InitAFC();
BrAFC

%% Creating a Falsification Problem
% Given a requirement R and some parameter range, we want to find a
% parameter value for which the system violates R.

%% 
% First we create the parameter search domain and load specifications. 
AFC_Falsify = BrAFC.copy();
AFC_Falsify.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 20; 10 60]);
STL_ReadFile('AFC_simple_spec.stl');
req = BreachRequirement(AF_alw_ok);

%% 
% Then we create the falsification problem and solve it.
falsif_pb = FalsificationProblem(AFC_Falsify, req);
falsif_pb.solve();

%% 
% That's it. The default solver found a trace violating the specification
% AF_alw_ok. 

%% Getting and Plotting the Falsifying Trace
AFC_False = falsif_pb.GetFalse(); % AFC_False contains the falsifying trace
AFC_False.PlotSignals({'Pedal_Angle','Engine_Speed','AF'});

%% Some figure cosmetics
subplot(3,1,3); set(gca, 'XLim', [10 40]);
plot([0 41], (1+0.005)*[14.7 14.7],'r');
plot([0 41], (1-0.005)*[14.7 14.7],'r');


%% Getting and Plotting the Falsifying Trace (ct'd)
% We can also examine more closely at the violation by plotting the
% robust satisfaction function for the formula and its subformulas:
AFC_False.PlotDiagnostics(); 

%% Trying a Harder Falsification Instance
% We were able to falsify AF_alw_ok which says that AF should stay within 1% of
% AFref. What if we soften the requirement and make it 5% ? 

req.SetParam('tol', 0.05);
falsif_pb2 = FalsificationProblem(AFC_Falsify, req);
res = falsif_pb2.solve();

%%
% The solver failed, we might want to understand what happened. 

%% Understanding the Solver Behavior
% The solver has reached the maximum number of objective function
% evaluation (i.e., number of simulations and property evaluation in our
% case). We can allocate more.  
falsif_pb2.max_obj_eval = 1000; % give it 1000 evalulations - default is 300

%%
% The solver can also have a maximum time to run (default is inf). 
falsif_pb.max_time


%% Understanding the Solver Behavior (ct'd)  
% The default solver follows a strategy alternating different phases. 
% # Corners: computes the objective function for corner values  
% # Global (quasi-random): computes the objective function for a number of trials in the
% parameter domain using quasi-random sampling  
% # Local: sorts the results obtained in the global phase and perform
% local optimizations using Nelder-Mead algorithm starting with the best
% value of the last quasi-random phase
%

%% Understanding the Solver Behavior (ct'd) 
% We can examine which parameter value was tested as follows: 

BrLog = falsif_pb2.GetLog(); 
BreachSamplesPlot(BrLog);

%% Understanding the Solver Behavior (ct'd) TODO META
% To control the default solver behavior, the following options can be set: 
falsif_pb2.solver_options

%% 
% * num_corners is the number of corners tested during the corners phase
% * num_quasi_random_samples is the number of quasi-random samples tested
% during the global phase
% * local_max_obj_eval is the maximum number of objective evaluations
% during one local phase
% * local_options allows to fine-tune the local optimizer. It can be
% changed using the Matlab optimset function


%% Understanding the Solver Behavior (ct'd)
% Falsification can be resumed by running solve() again. The main advantage
% is that it will skip already tried trials in the global phase. In particular,
% corner values won't be tried again. 
falsif_pb2.solve(); % resume falsification - will run for 1000 simulations

%% Understanding the Solver Behavior (ct'd) 
% We can check that more values have been explored:  
BrLog2 = falsif_pb2.GetLog(); 
BreachSamplesPlot(BrLog2);

%% Other Solvers 
% We can get a list of available solvers by
BreachProblem.list_solvers() 

%%
% Each of these solvers can be setup with default options using 
% setup_solver(). Some require the Optimization toolbox.

%%
% E.g, CMA-ES is a popular evolutionary algorithm included in the list of
% solvers.
 
falsif_pb3 = FalsificationProblem(AFC_Falsify, req,{'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 20; 10 60]);
falsif_pb3.setup_solver('cmaes');

%% Running CMAES
falsif_pb3.solve();

%% Customizing/Interfacing a New Solver
% New solvers can be interfaced (or new interface created) by creating a new class:  
type MyFalsificationProblem.m

%% Customizing/Interfacing a New Solver (ct'd) 
myfalsif_problem = MyFalsificationProblem(AFC_Falsify, req);
myfalsif_problem.solve();







