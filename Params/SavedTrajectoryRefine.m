function P = SavedTrajectoryRefine(P, paramValues)
% This is a "Refine" method meant to load correct parameter variables when
% a trajectory has been saved.

nb = 1; % Number of samples, should be 1!
P.epsi = repmat(P.epsi,[1 nb]);

P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = paramValues;

% manage traj_ref and traj_to_compute: we make the supposition that no
% generated parameter vector is equal to an previously computed parameter
% vector
[~,P.traj_to_compute] = unique(P.pts(1:P.DimP,:)','rows','first');
P.traj_ref = zeros(1,size(P.pts,2));
P.traj_to_compute = sort(reshape(P.traj_to_compute,1,[]));

end

