function P = HaltonRefine(P, nb, varargin)
%HALTONREFINE samples quasi-uniformly a parameter set using Halton
% sequence. This function should be used directly as it does not manage the
% fields traj_ref and traj_to_compute, but rather through QuasiRefine
% function.
%
% Synopsis:  P = HaltonRefine(P, nb[, step][, 'strictlyInside'])
%
% Inputs:
%  - P                : a parameter set. It may contain many parameter
%                       vectors
%  - nb               : the number of parameter vectors generated for each
%                       initial parameter vector
%  - step             : (optional, default=1)
%  - 'strictlyInside' : (optional, not set by default) If the string
%                       'striclyInside' is provided, the generated
%                       parameter sets are such that all new boxes are
%                       stricly inside the initial one. Otherwise, the
%                       centers of the new boxes are in the initial one,
%                       but some part of the boxes may overtake the initial
%                       one.
%
% Output:
%  - P : the new parameter set
%
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys); % Create default parameter set for system Sys
%   Ph = QuasiRefine(P, 1000,'halton'); % Sample with 1000 points
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph);   % plots the generated points
%
% Credit:  John Burkardt, 2003
%
%See also QuasiRefine Refine RandomLogRefine LogNRefine
%

if(nb<=1)
    return;
end

if(nargin==2)
    step = 1;
    strictlyInside = false;
elseif(nargin==3)
    if ischar(varargin{1})
        strictlyInside = strcmpi(varargin{1},'strictlyinside');
        step = 1;
    else
        strictlyInside = false;
        step = varargin{1};
    end
else
    step = varargin{1};
    strictlyInside = strcmpi(varargin{2},'strictlyinside');
end


dim_num = numel(P.dim);

seed = 1*ones(dim_num,1);
leap = ones(dim_num,1);

base = primes(dim_num*dim_num+1);
base = base(1:dim_num);

r = halton_sequence(dim_num,nb, step, seed, leap, base);
r = kron(r, ones(1,size(P.pts,2)));


old_epsi = P.epsi;
new_epsi = P.epsi/(nb^(1/dim_num));
P.epsi = repmat(new_epsi,[1 nb]);
% old version
%P.epsi = repmat(P.epsi,[1 nb])/(floor(nb^(1/dim_num)));

if(strictlyInside)
    width = 2*(old_epsi - new_epsi);
    mini = P.pts(P.dim,:) - (old_epsi - new_epsi);
else
    width = 2 * old_epsi;
    mini = P.pts(P.dim,:) - old_epsi;
end
P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(mini,[1 nb]);

if isfield(P,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

P = Preset_traj_ref(P);

end

function r = halton_sequence ( dim_num, n, step, seed, leap, base )
% was originally I4_TO_HALTON_SEQUENCE: N elements of an DIM_NUM-dimensional Halton sequence.
%
%  Author:John Burkardt
%

dim_num = floor ( dim_num );
n = floor ( n );
step = floor ( step );
seed(1:dim_num) = floor ( seed(1:dim_num) );
leap(1:dim_num) = floor ( leap(1:dim_num) );
base(1:dim_num) = floor ( base(1:dim_num) );
r(1:dim_num,1:n) = 0.0;

for i = 1: dim_num
    
    seed2(1:n) = seed(i) + step * leap(i) : leap(i) : ...
        seed(i) + ( step + n - 1 ) * leap(i);
    
    base_inv = 1.0 / base(i);
    
    while ( any ( seed2 ~= 0 ) )
        digit(1:n) = mod ( seed2(1:n), base(i) );
        r(i,1:n) = r(i,1:n) + digit(1:n) * base_inv;
        base_inv = base_inv / base(i);
        seed2(1:n) = floor ( seed2(1:n) / base(i) );
    end
end

end
