function S = SXf2X0(S0)
  
  DimX = S0.DimX;
  S = S0;
  S.pts(1:DimX,:) = S.Xf;
  S = SPurge(S);
  
  
  