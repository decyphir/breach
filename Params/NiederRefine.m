function Sh = NiederRefine(S, nb, step)
% NIEDERREFINE  Sample quasi-uniformly a parameter set using niederreiter2 sequence
% 
% Synopsis:  Ph = NiederRefine(S, nb)
%  
% Example: 
%
%   CreateSystem;
%   P = CreateParamSet(Sys); % Create default parameter set for system Sys
%   Ph = NiederRefine(P, 1000); % Sample with 1000 points
%  
%   SplotBoxPts(P); % Parameter set before sampling 
%   SplotPts(Ph);   % plots the generated points
%
% Credit:  John Burkardt, 2003
%  
  
  dim_num = numel(S.dim);
  
  if (nargin==2)
    step=1;
  end
  
  seed = 1*ones(dim_num,1);
  leap = ones(dim_num,1);
  
  base = primes(dim_num*dim_num+1)
  base = base(1:dim_num);
  
  r = niederreiter2_generate ( dim_num, nb, seed )
  
  r = kron(r, ones(1,size(S.pts,2)));
  
  A = 2*S.epsi;
  a = S.pts(S.dim,:)-S.epsi;
  
  Sh = S;
  Sh.pts = repmat(S.pts,[1 nb]);
  
  Sh.pts(S.dim,:) = repmat(A,[1 nb]).*r+repmat(a,[1 nb]);

  %Sh.epsi = repmat(Sh.epsi,[1 nb])/(floor(nb^(1/dim_num)));

  Sh.epsi = kron(Sh.epsi,ones(1,size(Sh.pts,2)))/(floor(nb^(1/dim_num)));

  %Sh.epsi = kron(Sh.epsi, ones(1,size(Sh.pts,2)))/nb;
  
  if (isfield(S,'selected'))
    Sh.selected = zeros(1, size(S.pts,2));    
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
