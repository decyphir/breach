function val = DistFromRange(traj, SignalRanges, t)

if exist('t', 'var')
    Xt = interp1(traj.time', traj.X', t')';
else
    Xt = traj.X;
end

minXt = SignalRanges(:,1);
maxXt = SignalRanges(:,2);

maxXt = repmat(maxXt, 1, size(Xt,2)); 
minXt = repmat(minXt, 1, size(Xt,2)); 

dist_maxX = maxXt-Xt;
dist_minX = Xt-minXt;

% takes the minimum, component-wise
dist_signals = ((dist_maxX + dist_minX) - abs(dist_maxX - dist_minX))/2;
val = min(dist_signals, [], 1);

end






















