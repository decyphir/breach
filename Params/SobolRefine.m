function Sh = SobolRefine(S, nb, step)
% SOBOLREFINE  Sample quasi-uniformly a parameter set using Sobol sequence
% 
% Synopsis:  Ph = SobolRefine(S, nb)
%  
% Example: 
%
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = SobolRefine(P, 1000); % Sample with 1000 points
%  
%   SplotBoxPts(P); % Parameter set before sampling 
%   SplotPts(Ph);   % plots the generated points
%
% Credit:  John Burkardt, 2006
%  
  
  dim_num = numel(S.dim);
  
  if (nargin==2)
    step=0;
  end
      
  r = i4_sobol_generate ( dim_num, nb, step );
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
  