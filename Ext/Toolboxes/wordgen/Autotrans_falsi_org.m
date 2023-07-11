%% STILL UNDER DEBUG

%% This script set up the falsification problem 
%% for the Autotransmission model

close all
clear all
warning('off', 'ALL')

thao =0; 
if thao
    addpath('/Users/thaodang/Metaheuristics/src')
    addpath('/Users/thaodang/Metaheuristics/breach-dev')
    addpath('.')    
    InitBreach('/Users/thaodang/Metaheuristics/breach-dev',true); % forces initialization from folder in Metaheuristics
else
    addpath('../../Metaheuristics/src')
end
model_name = 'autotrans_mod04';
fprintf('\n Creating breach interface with simulink model %s\n',model_name)

simTime = 30 ; 
fprintf('\n Simulation time horizon is %d seconds\n',simTime)

IO_signal_names = {'In1','In2','RPM','gear','speed'};

%% Specifying STL formula
%fprintf('\n The STL formula is\n ')
f1 = STL_Formula('f1','alw(RPM[t]<2520)');
f2 = STL_Formula('f2','ev_[0,10](speed[t]>50)');
%phi = STL_Formula('phi','alw_[0,10](RPM[t]<2520)');
phi = STL_Formula('phi','not(f1 and f2)');
%phi = STL_Formula('phi','alw_[0,10](speed[t]<50)');

% Can rewrite formula: [f1 => phitest2] where
% phitest1 = STL_Formula('phitest1','ev(RPM[t]>2520)'); % not(f1)
% phitest2 = STL_Formula('phitest2','alw_[0,10](speed[t]<50)'); % not(f2)


%% Set input signals

fprintf('\n Parametrizing input signals throttle and break....\n')
fprintf('\n Input signals are parametrized as piecewise constant\n ')

nb_ctr_pts = [ 7; 3 ];
input_ranges = [ 35 100; 0 40 ]; 
%signal_types = { 'UniStep', 'UniStep' };
%%%% thao signal_types = { 'VarStep', 'VarStep' };
input_signal_names = {'In1','In2'};
signal_gen_method = {'previous','previous'}; 


%%%%
% Insig = load(sigfilename, '-ascii');
% scaling=0.65e-7; 
%     %%Passband 0.4e-7 Oui les deux -- 0.5e-7 Oui les deux -- 
%     %0.6e-7 discrepance 19 Oui, uni NON -- %0.8e-7 NON les deux
%     
%     %%Lowpass
%     %%0.4e-4;
%     
% Insig(:,1) = scaling*In1(:,1);
% time = Insig(:,1);
%%%% signal_gen_method = {'linear','linear'}; 
%%%% In1_u0 = 
%%%% In1_dt0 = time(2) - time(1)
%%%% pas de dernier interval dt 


%from old code
%input_gen.type = 'UniStep';
%input_gen.method = {'previous','previous'};

% grid size collumn on the range of each input signal
gridsize_vector = [ 4; 4 ];


%%%% Once the above system specifications and falsification options are given,
%%%% the following part of the code need not be modified by the user
MetaObj = MetaFalsify(model_name,IO_signal_names);
% specify the simulation time
MetaObj.SimTimeSetUp(simTime);

% specify the class of input signals
MetaObj.InputSignalSetUp(input_signal_names,signal_gen_method,nb_ctr_pts,input_ranges);

% specify the property to falsify
MetaObj.STLFormulaSetUp(phi);

% set up a grid on the input ranges, to estimate coverage
MetaObj.GridSetUp(gridsize_vector,nb_ctr_pts);        


%% Start the falsification process
%%%[r,falsified,total_nb_sim,falsi_point] = MetaObj.MetaCall();

%%% Set up falsification options

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
    MetaObj.re_init_strategy = 2; %2; 
    % re_init_strategy=0 to pick randomly from the whole space
    % re_init_strategy=1 to pick randomly from xlog
    % re_init_strategy=2 to pick randomly from xbest
            
    % re_init_num_xbest: window of choice from xbest, for picking initial point         
    MetaObj.re_init_num_xbest = 200;

    % num_solvers=nb of solvers %%other than pseudorandom sampling
    % TODO add solver_list, and init num_solver as numel(solver_list)
    MetaObj.num_solvers=4; 
    
    %% limit on nb of solver calls
    MetaObj.nb_solver_calls = 2  %1 %30
    
    MetaObj.start_solver_index = 2; %3; %1; %PR 0, cmaes 1, SA 2, GNM 3 
    
    MetaObj.solver_time = [500 2000 500 900];
    MetaObj.max_obj_eval = [ 1000 20000 1000 50000 ];
    MetaObj.seed = 5000;
    
    
fprintf('\n The falsification problem by metaheuristics is\n ')
MetaObj

%%%% Run the falsification
MetaObj.MetaCall();

%MetaObj.MetaSetupRun(Sys, phi)
%MetaObj.MetaSetupRun(MetaObj.Br, phi)

