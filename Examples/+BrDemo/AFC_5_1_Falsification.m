%% Falsification Problems

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
AFC_Falsify.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 20; 10 60]);
STL_ReadFile('AFC_simple_spec.stl');

%% 
% Then we create the falsification problem and solve it.
falsif_pb = FalsificationProblem(AFC_Falsify, AF_alw_ok);
falsif_pb.solve();

%% 
% That's it. The default solver found a trace violating the specification
% AF_alw_ok. 

%% Getting and Plotting the Falsifying Trace

AFC_False = falsif_pb.GetBrSet_False(); % AFC_False contains the falsifying trace
AFC_False.PlotSignals({'Pedal_Angle','Engine_Speed','AF'});

%% Some figure cosmetics
subplot(3,1,3); set(gca, 'XLim', [10 40]);
plot([0 41], (1+0.005)*[14.7 14.7],'r');
plot([0 41], (1-0.005)*[14.7 14.7],'r');


%% Getting and Plotting the Falsifying Trace (ct'd)
% We can also examine more closely at the violation by plotting the
% robust satisfaction function for the formula and its subformulas:
figure;AFC_False.PlotRobustSat(AF_alw_ok,3); % second argument is depth of formula decomposition

%% Trying a Harder Falsification Instance
% We were able to falsify AF_alw_ok which says that AF should stay within 1% of
% AFref. What if we soften the requirement and make it 5% ? 

AF_alw_ok2 = set_params(AF_alw_ok, 'tol', 0.05);
falsif_pb2 = FalsificationProblem(AFC_Falsify, AF_alw_ok2);
res = falsif_pb2.solve();

%%
% The solver failed, we might want to understand what happened. 

%% Understanding the Solver Behavior
% The solver stopped because a timeout was reached. We can allocate more
% time by changing max_time value:
falsif_pb2.max_time = 120; % give it two minutes - default is 60s

%%
% Note that the solver may also have a maximum number of objective function
% evaluation (i.e., number of simulations and property evaluation in our
% case).  
falsif_pb.max_obj_eval

%%
% Inf means there is no bound here. 

%% Understanding the Solver Behavior (ct'd)  
% The default solver follows a strategy with two phases:  
%
% # Global: computes the objective function for a number of trials in the
% parameter domain, starting with corners then quasi-random sampling  
% # Local:  sorts the results obtained in the global phase and perform
% local optimizations using Nelder-Mead algorithm starting with the best values as initial guesses  


%% Understanding the Solver Behavior (ct'd) 
% We can examine which parameter value was tested as follows: 

BrLog = falsif_pb2.GetLog(); 
BreachSamplesPlot(BrLog);

%% Understanding the Solver Behavior (ct'd)
% To control the default solver behavior, the following options can be set: 
falsif_pb2.solver_options

%% 
% * start_at_trial is the number of trials to skip. Can be considered
% a seed for random exploration, except that low numbers corresponds
% to corners of the parameter domain. It is 24 here because this is the
%    number of trials that were performed in the previous run. 
% * nb_new_trials is the number of trials to run in the next call to solve(). 
% * nb_local_iter is the number of Nelder-Mead iterations to perform in the
% local phase for each trial. 
% * local_optim_options allows to fine-tune the local optimizer. It can be
% changed using the Matlab optimset function


%% Understanding the Solver Behavior (ct'd)
% Falsification can be resumed by running solve() again. The main advantage
% is to skip already tried trials in the global phase. In particular,
% corner values won't be tried again. 
falsif_pb2.solve(); % resume falsification - will run at most 2 minutes

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
 
falsif_pb3 = FalsificationProblem(AFC_Falsify, AF_alw_ok2,{'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 20; 10 60]);
falsif_pb3.setup_solver('cmaes');

%% Running CMAES
falsif_pb3.solver_options.MaxFunEvals = 100; % max number of simulations 
falsif_pb3.solve();

%% Customizing/Interfacing a New Solver
% New solvers can be interfaced (or new interface created) by creating a new class:  
type MyFalsificationProblem.m

%% Customizing/Interfacing a New Solver (ct'd) 
myfalsif_problem = MyFalsificationProblem(AFC_Falsify, AF_alw_ok);
myfalsif_problem.solve();







