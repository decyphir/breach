%% This script set up the falsification problem 
%% for the Autotransmission model

close all
clear all
%warning('off', 'ALL')

addpath('/Users/thaodang/Metaheuristics/supp_code')
addpath('/Users/thaodang/Metaheuristics/src')
addpath('/Users/thaodang/Metaheuristics/breach-dev')
addpath('.')
        
InitBreach('/Users/thaodang/Metaheuristics/breach-dev',true); % forces initialization from folder in Metaheuristics

model_name = 'autotrans_mod04';
fprintf('\n Creating breach interface with simulink model %s\n',model_name)

simTime = 30 ; 
fprintf('\n Simulation time horizon is %d seconds\n',simTime)

IO_signal_names = {'In1','In2','RMP','gear','speed'};
BrSys = CoverageBreachSet(model_name,{},[],IO_signal_names);
%BrSys = CoverageBreachSet(model_name,{});
BrSys.SetTime([0 simTime]);

%% Set input signals

fprintf('\n Parametrizing input signals throttle and break....\n')
fprintf('\n Input signals are parametrized as piecewise constant\n ')
input_gen.type = 'UniStep';

N1 = 7; N2 = 3;
fprintf('Number of control points for throttle input is %d\n',N1)
fprintf('Number of control points for break input is %d\n',N2 )
input_gen.cp = [N1 N2];
input_gen.method = {'previous','previous'};
BrSys.SetInputGen(input_gen);

%% Specifying parameter names
for i=0:N1-1
    signal_u0{1,i+1}=strcat('In1_u',num2str(i));
end

for i=0:N2-1
     signal_u1{1,i+1}=strcat('In2_u',num2str(i));    
end
signal = [signal_u0,signal_u1];

%% Initializing CBS object parameters

% Input ranges
% fprintf('\n Range of throttle is [35,100]\n')
% fprintf('Range of break is [0,40] \n')
% fprintf('\n Grid discretization unit for both signal ranges is 4 units\n')
R1 = [35,100];
R2 = [0, 40];
Sys = BrSys.copy();
%signal
Sys.SetParamRanges(signal,[ones(N1,1)*R1;ones(N2,1)*R2]);
Sys.SetEpsGridsize([4*ones(N1,1);4*ones(N2,1)]);
Sys.SetDeltaGridsize(2*Sys.epsgridsize);

%% Specifying STL formula
%fprintf('\n The STL formula is\n ')
f1 = STL_Formula('f1','alw(RPM[t]<2520)');
f2 = STL_Formula('f2','ev_[0,10](speed[t]>50)');
phi = STL_Formula('phi','not(f1 and f2)');
%phi = STL_Formula('phi','alw_[0,10](speed[t]<50)');

% Can rewrite formula: [f1 => phitest2] where
% phitest1 = STL_Formula('phitest1','ev(RPM[t]>2520)'); % not(f1)
% phitest2 = STL_Formula('phitest2','alw_[0,10](speed[t]<50)'); % not(f2)



MetaObj = MetaFalsify();




%%% Set up falsification options

    %% limit on nb of solver calls
    %this.nb_solver_calls = input('Specify Max Nb of Solver Calls: '); 
    MetaObj.nb_solver_calls = 30 %30 %1 %30
    %fprintf('\n Max Nb of Solver Calls: ',this.nb_solver_calls)

    %%%% Search Monitoring Parameters 
    %% cov_epsilon = input('Specify coverage increase threshold : '); 
    MetaObj.cov_epsilon = 1e-3;
    %% min robustness decrease in percentage
    MetaObj.rob_epsilon_percent = 0.05;
    %% min robustness stagnant monitoring window
    MetaObj.rob_stagnant_win = 1 
    %% coverage stagnant monitoring window
    MetaObj.cov_monitoring_win = 1;

    %%% Options for picking initial conditions
    MetaObj.re_init_strategy = 1; 
    % re_init_strategy=0 to pick randomly from the whole space
    % re_init_strategy=1 to pick randomly from xlog
    % re_init_strategy=2 to pick randomly from xbest
            
    % re_init_num_xbest: window of choice from xbest, for picking initial point         
    MetaObj.re_init_num_xbest = 1;

    % num_solvers=nb of solvers %%other than pseudorandom sampling
    % TODO add solver_list, and init num_solver as numel(solver_list)
    MetaObj.num_solvers=4; 
    
    MetaObj.start_solver_index = 3; %1; %PR 0, cmaes 1, SA 2, GNM 3 
    
    MetaObj.solver_time = [200 1200 200 200];
    MetaObj.max_obj_eval = [ 2000 2500 2000 2000 ];
    MetaObj.seed = 5000;
    
    
fprintf('\n The falsification problem by metaheuristics is\n ')
Sys

%%%% Run the falsification
MetaObj.MetaSetupRun(Sys, phi);


