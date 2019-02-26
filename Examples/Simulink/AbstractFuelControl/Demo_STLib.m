close all;
clear;

%%
STL_ReadFile('stlib.stl');

%%
SG = SignalGenerator()
%%
% Mixed bag of signals
% Here we generate 1 random signals for each signal configuration 
SG.GenRandomSignalsFromCfg(1:26, 5); 
SG.PlotSignals
%%
%phi = req_stable
%phi = req_inc
%phi = req_dec
%phi = req_ev_ring_once
phi = ev_ring_twice
%phi = req_overshoot

% checks the stable template on all those signals
val = SG.GetRobustSat(phi);

% reveals positive and negative 
idx_pos = find(val>0);
idx_neg = find(val<0);

%%
if ~isempty(idx_neg)
    figure;
    SG.PlotSignals([],idx_neg )
    set(gcf, 'Name','Negative');
end

%%
if ~isempty(idx_pos)
    figure;
    SG.PlotSignals([],idx_pos )
    set(gcf, 'Name','Positive');
end

%% 
configs = GetParam(SG.P, 'config');
pos_types = {};
for idx= idx_pos
    icfg  = configs(idx);    
    pos_types = {pos_types{:}, SG.configs{icfg}.name}; 
end
pos_types = unique(pos_types)'

configs = GetParam(SG.P, 'config');
neg_types = {};
for idx= idx_neg
    icfg  = configs(idx);    
    neg_types = {neg_types{:}, SG.configs{icfg}.name}; 
end
neg_types = unique(neg_types)'







