function X =  traj2pts(S)
%
% TRAJ2PTS converts trajectories to set of points
%

  X = cat(2,S.traj{1}.X);
