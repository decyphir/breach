function S2 = SCopyEmpty(S1)
% SCOPYEMPTY copies a parameter set, removing everything except dimension and  param names 
%
% Synopsis: Pnull = SCopyEmpty(P)
%
  
  S2.DimX = S1.DimX;
  S2.DimP = S1.DimP;
  S2.dim  = S1.dim;
  S2.ParamList = S1.ParamList;
  S2.pts = [];
  S2.epsi = [];
  
  