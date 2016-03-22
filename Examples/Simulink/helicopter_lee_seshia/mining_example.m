%% Analysis of a simple Simulink example (Example 2.5 from Lee and Seshia, Chapter 2)
%
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

%%
% We extract the nominal parameter set and simulate the system using Breach
% command line interface.
Br.Sim();

%% 
% We plot the nominal trace using PlotSignals:
Br.PlotSignals();

%% Parameter synthesis
% In the following, we analyse how fast $\dot{\theta}$ converges toward 
% the input $\psi$ by using the parameter synthesis algorithm on a PSTL 
% formula. 

%% 
% The following loads the PSTL property phi in the workspace:
STL_ReadFile('specs.stl');

%% 
% We define the options for the parameter synthesis algorithm. The algorithm 
% works by performing a binary search over the possible values of the
% parameters.
prop_params.names = {'tau'};
prop_params.ranges = [0 10];  

synth_pb = ParamSynthProblem(Br, phi, prop_params.names, prop_params.ranges); 
synth_pb.solve();

%%
% The value computed by solving the synthesis problem is store in x_best:
tau = synth_pb.x_best;

%% 
% We update the formula and plot its satisfaction:
phi_tight = set_params(phi, 'tau', tau);
Br.PlotRobustSat(phi_tight)

%% Falsification 
% Now we try to find a configuration of the system which will violate the 
% previous specification, satisfied by the nominal trajectory.

%% 
% We defines the system parameter name(s) and range(s).  
falsif_params.names = {'K' ... ,
                      };
%% 
% The parameter K can vary between 9 and 11
falsif_params.ranges = [9 11; ...
];

%% 
% We prepare a falsification problem.
falsif_pb = FalsificationProblem(Br, phi_tight, falsif_params.names, falsif_params.ranges)
falsif_pb.solve();

%%
% We plot the falsifying trajectory.  
BrFalse = falsif_pb.GetBrSet_False();
BrFalse.PlotRobustSat(phi_tight)


%% Mining  
% Next, we use the ReqMining problemroutine to combine parameter synthesis and
% falsification until finding an STL formula that the falsifier cannot
% violate. 
mining_phi = ReqMiningProblem(Br,phi, falsif_params, prop_params);
mining_phi.solve();


%% Validation 
% The model analysed in this example is quite trivial, and we know that the
% longest time to converge toward the input is obtained by the smallest
% value for K. Hence we check that the value found before is conservative,
% i.e., is larger than the one obtained with the smallest value of K.






