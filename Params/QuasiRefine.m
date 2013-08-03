function Ph = QuasiRefine(P, nb, varargin)
%QUASIREFINE Sample quasi-uniformly a parameter set. The number of
% uncertain parameters must be lower or equal to 40 (when using Sobol).
%
% Synopsis:  Ph = QuasiRefine(P, nb [,step] [, 'strictlyInside'] [, algo] )
%
% Inputs:
%  - P    : The parameter set to refine. May contains many set of parameter
%           values
%  - nb   : how many parameter set are generated for each set of parameter
%           values. If lower or equal to one, nothing is done. If nb is not
%           an integer, it is rounded toward 0.
%  - step : number of generated parameter to skip (optional, default=0)
%  - strictlyInside : (Optional, default=not set) If set (aka, write
%                     'strictlyInside' in the function call), all generated
%                     parameter sets are strictly inside the initial one.
%                     Otherwise, the center of each generated parameter set
%                     are inside the initial one, but the associated boxes
%                     may overtake the initial one. TODO: THE
%                     strictlyInside OPTION IS NOT AVAILABLE WITH THE SOBOL
%                     ALGORITHM.
%  - algo : (Optional, default='Halton') String defining which algorithm
%           should be used to generate random number used to generate new
%           parameter sets. Can be either 'Halton' or 'Sobol'.
%
% Output:
%  - Ph : The new parameter set
%
% Example:
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = QuasiRefine(P, 1000); % Sample with 1000 points
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph);   % plots the generated points
%
%See also RandomLogRefine LogNRefine Refine
%

% we manage the inputs
if(nb<=1)
    return;
end
nb = floor(nb);

if(nargin==2) % no optional parameter
    step = 1;
    strictlyInside = false;
    algo = 'halton';
elseif(nargin==3) % one optional parameter
    if ischar(varargin{1})
        strictlyInside = strcmpi(varargin{1},'strictlyinside');
        if isMethodValid(varargin{1})
            algo = varargin{1};
        else
            algo = 'halton';
        end
        step = 1;
    else
        strictlyInside = false;
        algo = 'halton';
        step = varargin{1};
    end
elseif(nargin==4) % two optional parameters
    if ischar(varargin{1})
        step = 1;
        strictlyInside = strcmpi(varargin{1},'strictlyinside');
        if isMethodValid(varargin{2})
            algo = varargin{2};
        else
            algo = 'halton';
        end
    else
        step = varargin{1};
        strictlyInside = strcmpi(varargin{2},'strictlyinside');
        if isMethodValid(varargin{2})
            algo = varargin{2};
        else
            algo = 'halton';
        end
    end
else % all optional parameters
    step = varargin{1};
    strictlyInside = strcmpi(varargin{2},'strictlyinside');
    if isMethodValid(varargin{3})
        algo = varargin{3};
    else
        algo = 'halton';
    end
end

if(strcmpi(algo,'sobol') && numel(P.dim)>40)
    warning('QuasiRefine:InappropriateAlgo',...
        'The sobol algorithm is usable up to dimension 40, switched to halton algorithm.');
    algo = 'halton';
end

if strcmpi(algo,'sobol')
    if (strictlyInside)
        Ph = SobolRefine(P,nb,step,'strictlyInside');
    else
        Ph = SobolRefine(P,nb,step);
    end
else % default = halton
    if (strictlyInside)
        Ph = HaltonRefine(P,nb,step,'strictlyInside');
    else
        Ph = HaltonRefine(P,nb,step);
    end
end

X = Ph.pts(1:Ph.DimP,:)';
[~,IA,IC] = unique(X,'rows');

Ph.traj_ref = IC';
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

function valid = isMethodValid(str)
valid = ( strcmpi(str,'halton') || strcmpi(str,'sobol') );
end

