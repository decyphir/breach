function P = SPurge(P)
%SPURGE removes all fields related to a specific computation of
% trajectories. This function does not remove properties values.
% 
% Synopsis: P = SPurge(P)
% 
% Input:
%  - P : the parameter set to purge
% 
% Output:
%  - P : the parameter set after cleaning
% 
%See also SPurge_props ComputeTraj
%

if isfield(P,'traj')
    P = rmfield(P, 'traj');
end
if isfield(P,'Xf')
    P = rmfield(P, 'Xf');
end

if isfield(P,'XSf')
    P = rmfield(P, 'XSf');
end
if isfield(P,'XS0')
    P = rmfield(P, 'XS0');
end
if isfield(P,'ExpaMax')
    P = rmfield(P, 'ExpaMax');
end

% Reset field traj_to_compute and traj_ref
[~,P.traj_to_compute] = unique(P.pts(1:P.DimP,:)','rows');
P.traj_to_compute = sort(reshape(P.traj_to_compute,1,[]));
P.traj_ref = zeros(1,size(P.pts,2));

end
