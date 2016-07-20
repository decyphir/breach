%% Analysis of a simple Simulink example (Example 2.5 from Lee and Seshia, Chapter 2)
% 

%% Interfacing and computing a nominal trajectory
% First, we interface a model with Breach and use Breach interface to
% compute a nominal simulation. 

%%
% Breach initialization. Make sure Breach home is in your Matlab path.  
clear 
InitBreach;

%%
% We initialize the system interface with Breach. This creates an object  
% called Br representing this interface. Note that Breach also creates 
% its own copy of the Simulink model. 
init_helicopter;
Br.PrintAll;

%%
% Simulate the system.
Br.Sim();

%% 
% We plot the nominal trace using PlotSignals:
Br.PlotSignals();

%% Parameter synthesis
% In the following, we analyse how fast $\dot{\theta}$ converges toward 
% the input $\psi$ by using the parameter synthesis algorithm on a PSTL 
% formula. 

%% 
% The following loads the PSTL property phi in the workspace.
STL_ReadFile('specs.stl');

%%
% It's defined as: 
% param tau = 5 
% phi :=  ev_[0,tau] ( abs(theta_dot[t]-psi[t]) < 0.01 ) 
%
% We can check if it is satisfied and plot the satisfaction functions,
% Boolean and quantitative 
Br.CheckSpec(phi) 
Br.PlotRobustSat(phi)


%% 
% We define the options for the parameter synthesis algorithm. The algorithm 
% works by performing a binary search over the possible values of the
% parameter until finding the smallest value for which phi is satisfied.
prop_params.names = {'tau'};
prop_params.ranges = [0 10];  

synth_pb = ParamSynthProblem(Br, phi, prop_params.names, prop_params.ranges); 
synth_pb.solve();

%%
% The value computed by solving the synthesis problem is store in x_best:
tau_tight = synth_pb.x_best;

%% 
% We update the formula and plot its satisfaction:
phi_tight = set_params(phi, 'tau', tau_tight);
Br.PlotRobustSat(phi_tight)

%% Falsification 
% Now we try to find a configuration of the system which will violate the 
% previous specification, satisfied by the nominal trajectory.

%% 
% We defines the system parameter name(s) and range(s).  Here, we consider
% that the parameter K can vary between 9 and 11.
falsif_params.names = {'a' ... ,
                      };
falsif_params.ranges = [0.9 1.1; ...
];

%% 
% We prepare a falsification problem.
falsif_pb = FalsificationProblem(Br, phi_tight, falsif_params.names, falsif_params.ranges);
falsif_pb.solve();

%%
% We plot the falsifying trajectory.  
BrFalse = falsif_pb.GetBrSet_False();
BrFalse.PlotRobustSat(phi_tight)


%% Mining  
% Next, we use the ReqMiningProblem routine to combine parameter synthesis and
% falsification until finding an STL formula that the falsifier cannot
% violate. 
mining_phi = ReqMiningProblem(Br,phi, falsif_params, prop_params);
mining_phi.solve();


%% Validation 
% The model analysed in this example is quite trivial, and we know that the
% longest time to converge toward the input is obtained by the smallest
% value for K. This is what the solver found in the first iteration: K=9.
% Note that in the second iteration, the solver didn't retry K=9. 




