TA_filename = 'SimpleUppaalTA.xml';
num_evts = 10;
signals = {'sig1','sig2'};


labels_ranges.all__.sig1 = [5 20];   
labels_ranges.all__.sig2 = [-10 0];

labels_ranges.phase1.sig2 = -4; % sig2 fixed in phase1

labels_ranges.phase2.sig1 = 10; % sig1 constant in phase2

B_uppaal = BreachTASignalGen( ...
     signals, ...     % signal names
     TA_filename, ...  % automata file
     labels_ranges, ... % mapping between labels of time words and ranges/values for signals
     num_evts);          % length of words in number of labels or events

B_uppaal.SampleDomain(10);
B_uppaal.Sim(0:.01:50);

B_uppaal.PlotSignals()