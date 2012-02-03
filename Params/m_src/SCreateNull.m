function S = SCreateNull(S)
  
  if (nargin ==0)
    S=[];
  end
  
  S.pts = [];
  S.traj = [];
  S.epsi = [];
  S.ToCheck = [];
  S.XS0 = [];
  S.Xf = [];
  S.ExpaMax = [];
  S.XSf = [];    