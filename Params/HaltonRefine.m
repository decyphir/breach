function P = HaltonRefine(P, nb, step)
% HALTONREFINE  Sample quasi-uniformly a parameter set using Halton sequence
%
% Synopsis:  P = HaltonRefine(P, nb [, step] )
%
% Input:
%  - P    : 
%  - nb   : 
%  - step : 
%
% Output:
%  - P : 
%
% Example:
%
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = HaltonRefine(P, 1000); % Sample with 1000 points
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph);   % plots the generated points
%
% Credit:  John Burkardt, 2003
%

if(nb<=1)
    return;
end

dim_num = numel(P.dim);

if(nargin==2)
    step=1;
end

seed = 1*ones(dim_num,1);
leap = ones(dim_num,1);

base = primes(dim_num*dim_num+1);
base = base(1:dim_num);

r = halton_sequence(dim_num,nb, step, seed, leap, base);
r = kron(r, ones(1,size(P.pts,2)));

width = 2*P.epsi; % NM: TODO : update such that all new parameter set are included
mini = P.pts(P.dim,:)-P.epsi;  % into the initial ones

P.pts = repmat(P.pts,[1 nb]);

P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(mini,[1 nb]);

P.epsi = repmat(P.epsi,[1 nb])/(nb^(1/dim_num));
%P.epsi = repmat(P.epsi,[1 nb])/(floor(nb^(1/dim_num)));
%P.epsi = repmat(P.epsi, [1 nb])/nb;

if isfield(P,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

end

function r = halton_sequence ( dim_num, n, step, seed, leap, base )

%% was originally I4_TO_HALTON_SEQUENCE: N elements of an DIM_NUM-dimensional Halton sequence.
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
