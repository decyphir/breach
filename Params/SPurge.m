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

P = Preset_traj_ref(P);

end
