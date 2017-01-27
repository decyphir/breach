function S = SXf2X0(S0)
% SXF2X0 creates a new parameter set from the end points of trajectories
%
% Synopsis: Pnew = SXf2X0(Pf)  
%  
%  
 
  DimX = S0.DimX;
  S = S0;
  S.pts(1:DimX,:) = S.Xf;
  S = SPurge(S);
  S = SPurge_props(S);  
  
