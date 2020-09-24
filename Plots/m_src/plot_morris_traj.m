function plot_morris_traj(S)
%  
%  PLOT_MORRIS_TRAJ plots Morris trajectories 
%
%  Useful to illustrate the Morris method
%  
%  Ex:
%    Pt.dim = 1:3;
%    Pt.epsi = [.5 .5 .5]';
%    Pt.pts = [.5 .5 .5]';
%    Ptr = pRefine(Pt, 4, 3);
%    plot_morris_traj(Ptr);
%
  

  n = numel(S.dim);
  sz = size(S.pts,2);
  r = sz/(n+1);
  
  SplotPts(S,[],[],{'.r', 'MarkerSize',16}); 
  for i = 0:r-1
    SplotPts(S,[],i*(n+1)+1:(i+1)*(n+1),{'-', 'MarkerSize',16, 'LineWidth', 2});
  end
    
   
