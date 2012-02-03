function plot_morris_traj(S)
%  
%  PLOT_MORRIS_TRAJ plots Morris trajectories 
%
%  Useful to illustrate the Morris method
%  
%  Ex.: 
% 
%  S = SCreate(ones(3,1), ones(3,1))
%  Sr = pRefine(S, 10, 3)  % 10-grid  with 3 trajectories
%  

  n = numel(S.dim);
  sz = size(S.pts,2);
  r = sz/(n+1);
  
  SplotPts(S,[],[],{'.r', 'MarkerSize',16}); 
  for i = 0:r-1
    SplotPts(S,[],i*(n+1)+1:(i+1)*(n+1),{'-', 'MarkerSize',16});
  end
    
   