function S = VoronoiRefine(S0)
% VORONOIREFINE Voronoi diagram based refinement (experimental)
%  
%   P = VoronoiRefine(P0)
% 
%   Returns a sampling P refining P0 based on the Voronoi diagram of
%   P0.
    
  n = numel(S0.dim);
  
  S.dim = S0.dim;
  S.pts =[];
  S.epsi = [];
  
  opt = {'Qt', 'Qbb','Qc', 'Qz'}; 
  
  X = S0.pts(S0.dim,:)'; 
  Max = max(X)+max(S0.epsi');
  Min = min(X)-max(S0.epsi');
  
  [V C] = voronoin(X, opt );
  
  V = unique(V, 'rows');
  
  % keep only vertices inside the initial domain
  
  for ix = 1:size(X, 2)
    iin = find(V(:,ix)<=Max(ix));
    V = V(iin,:);
  
    iin = find(V(:,ix)>=Min(ix));
    V = V(iin,:);
  end
 
  S.pts = [S0.pts repmat(S0.pts(:,1), [1 size(V,1)-1])];
  nb0 = size(S0.pts,2);
  S.pts(S.dim, nb0+1:end) = V(2:end,:)';
  S.epsi = [S0.epsi/2 repmat(S0.epsi(:,1)/2, [1 size(V,1)-1])];  
  
  S.DimX= S0.DimX;
  S.DimP= S0.DimP;
  
