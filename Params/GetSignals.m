function X = GetSignals(traj,signals)
% GetSignals returns signal values for given name from a trace

if ischar(signals)
    signals = {signals};
end

[~, i_found, i_where] = intersect(signals, traj.signals, 'stable');
if numel(i_found) ~= numel(signals)
    i_not_found = setdiff(1:numel, i_found);
    error('Signal %s not found', signals{i_not_found(1)}); 
end

X = traj.X( i_where, : ); 


