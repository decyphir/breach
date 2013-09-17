function P = SPurge(P)
% SPURGE Remove all fields related to a specific computation of
% trajectories
%
%  Synopsis: P = SPurge(P)
%
%  Notes: does not remove properties values
%
%  SEE ALSO SPURGE_PROPS
%

try
    P = rmfield(P, 'traj');
end
try
    P = rmfield(P, 'Xf');
end

try
    P = rmfield(P, 'XSf');
end
try
    P = rmfield(P, 'XS0');
end
try
    P = rmfield(P, 'ExpaMax');
end

% Reset field traj_to_compute

X = P.pts(1:P.DimP,:)';
[~,IA,IC] = unique(X,'rows');

P.traj_ref = IC';
P.traj_to_compute = IA';

end
