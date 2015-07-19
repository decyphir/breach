%% Analysis of a simple Simulink example (Example 2.5 from Lee and Seshia, Chapter 2)
%

%% Interfacing and computing a nominal trajectory
% First, we interface a model with Breach and use Breach interface to
% compute a nominal simulation. 

%%
% Breach initialization. Make sure Breach home is in the path.  
clear 
InitBreach;

%%
% We initialize the system interface with Breach. This creates a structure 
% named Sys representing this interface. Note that Breach also creates 
% its own copy of the Simulink model. 
init_helicopter;

%%
% We extract the nominal parameter set and simulate the system using Breach
% command line interface.
P = CreateParamSet(Sys);
P = ComputeTraj(Sys,P); % P now contains a field P.traj with the result of 
                        % the simulation 

%%
% We plot the signals monitored by Breach:
figure;
SplotVar(P);

%% Parameter synthesis
% In the following, we analyse how fast $\dot{\theta}$ converges toward 
% the input $\psi$ by using the parameter synthesis algorithm on a PSTL 
% formula. 

%% 
% The following loads the PSTL property phi in the workspace:
QMITL_ReadFile('specs.stl');

%% 
% We define the options for the parameter synthesis algorithm. The algorithm 
% works by performing a binary search over the possible values of the
% parameters.
prop_opt.params = {'tau'};
prop_opt.monotony = [1];  % the formula is positively monotonic wrt tau 
prop_opt.ranges = [0 10];  
prop_opt.p_tol = [.01];   % precision of the binary search
                          

[p,~,P] = GetPropParamBin(Sys, phi, P, prop_opt);

%% 
% The following instantiate the value found by the parameter synthesis algorithm.
prop_opt.values= p;

%%
% We plot the satisfaction of the formula. 
figure;
SplotSat(Sys, P, phi, 3);

%% Falsification 
% Now we try to find a configuration of the system which will violate the 
% previous specification, satisfied by the nominal trajectory.

%% 
% We defines the system parameter name(s) and range(s).  
falsif_opt.params = {'K' ... ,
                    };
%% 
% The parameter K can vary between 9 and 11
falsif_opt.ranges = [9 11; ...
];

%% 
% We define the options for the falsifier optimization algorithm.  The 
% algorithm works by first sampling the parameter set with nb_init values, 
% then run a simulation for each of these values. If one is a falsifying simulation,
% it returns it, otherwise, it runs nb_iter iterations of Nelder Mead
% starting from each of them. It will stop if a falsifying trace is found
% or the algorithm performed nb_max_call simulations in total.
falsif_opt.nb_init = 10; 
falsif_opt.nb_iter = 10;
falsif_opt.nb_max_call = 100;

Pfalse = Falsify(Sys, phi, falsif_opt, prop_opt);

%%
% We plot the falsifying trajectory.  
figure;
SplotSat(Sys, Pfalse, phi, 3);

%% Mining  
% Next, we use the ReqMining routine to combine parameter synthesis and
% falsification until finding an STL formula that the falsifier cannot
% violate. 
iter_max= 10;
[p,~, Pr] = ReqMining(Sys, phi, falsif_opt, prop_opt, iter_max);
Psave(Sys, 'Pr'); % saves Pr in Breach GUI

%% Validation 
% The model analysed in this example is quite trivial, and we know that the
% longest time to converge toward the input is obtained by the smallest
% value for K. Hence we check that the value found before is conservative,
% i.e., is larger than the one obtained with the smallest value of K.
Pground = SetParam(P,'K', 9); 
Pground = ComputeTraj(Sys,Pground);
[p_ground,~, Pground] = GetPropParamBin(Sys, phi, Pground, prop_opt);




