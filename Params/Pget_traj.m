function traj = Pget_traj(P,i)


if isfield(P,'traj')&&numel(P.traj)>=i
    traj = P.traj{i};
else
    traj = {};
end

