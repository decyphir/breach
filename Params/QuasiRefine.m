function Ph = QuasiRefine(P, nb, step)
%QUASIREFINE  Sample quasi-uniformly a parameter set. The number of
% uncertain parameters must be lower or equal to 40 (when using Sobol).
%
% Synopsis:  Ph = QuasiRefine(P, nb [,step] )
%
% Input:
%  - P    : 
%  - nb   : 
%  - step : number of generated parameter to skip (default=0)
%
% Output:
%  - Ph :
%
% Example:
%
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = QuasiRefine(P, 1000); % Sample with 1000 points
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph);   % plots the generated points
%
% TODO optionalize the sequence used - so far Sobol
%
%See also RandomLogRefine LogNRefine Refine
%

if(nargin==2)
    step=0;
end

%if(numel(P.dim)<=40)
%    Ph = SobolRefine(P,nb,step);
%else
    Ph = HaltonRefine(P,nb,step);
%end

X = Ph.pts(1:Ph.DimP,:)';
[~,IA,IC] = unique(X,'rows');

Ph.traj_ref= IC';
Ph.traj_to_compute = IA';

if isfield(P,'traj')
    Ph.traj = P.traj;
    Ph.Xf = P.Xf;
    if ~isequal(P.pts(1:P.DimP,IA), vertcat(P.traj.param))
        Ph.traj_to_compute = IA';
    end
else
    Ph.traj_to_compute = IA';
end

end
