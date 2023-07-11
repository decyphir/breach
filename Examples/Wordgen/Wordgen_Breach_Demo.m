%% Import a Timed Automaton as a Breach Signal Generator
%  First we import a timed automaton defined with Wordgen web interface or prism, or Uppal.  
%  To configure the import, we need to specify the TA file, the desired
%  length of time words, and options for wordgen. Then we need to inform
%  Breach about how to convert labels into real values. 

TA_filename = 'driving_TA.prism';
num_evts = 10;
signals = {'throttle','brake'};

%% Define ranges for the signals for all labels
% The signals will take values in these ranges unless specified otherwise
% on a per label basis. 
labels_ranges.all__.throttle = [0 100];
labels_ranges.all__.brake = [0 350];


%% Definie specific ranges for each label, always starting with init label
% label init: Throttle takes default range but brake is 0  
labels_ranges.init.brake = 0;

% label a: no acceleration and no braking
labels_ranges.a.throttle =0; 
labels_ranges.a.brake =0; 

% label b: no acceleration, braking in range
labels_ranges.b.throttle =0; 

% label c: no braking
labels_ranges.c.brake =0; 

% label d: no acceleration
labels_ranges.d.throttle =0; 

% label e: no acceleration and no braking
labels_ranges.e.throttle =0; 

% label f: no acceleration and no braking
labels_ranges.f.throttle =0; 
labels_ranges.f.brake =0; 

% label g: no braking
labels_ranges.g.brake =0; 

% label h: no braking
labels_ranges.h.brake =0; 

%% Create BreachSignalGen object which can generate signals from the TA 
B_ta = BreachTASignalGen(signals, TA_filename,labels_ranges,num_evts);

%%  Uniform sampling
B_ta.SampleDomain(3)
B_ta.Sim(0:0.01:50);
B_ta.PlotSignals()




