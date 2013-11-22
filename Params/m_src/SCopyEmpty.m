function Pnull = SCopyEmpty(P)
% SCOPYEMPTY copies a parameter set, removing everything except dimension
% and parameter names
% 
% Synopsis: Pnull = SCopyEmpty(P)
%

Pnull.DimX = P.DimX;
Pnull.DimP = P.DimP;
Pnull.dim = P.dim;
Pnull.ParamList = P.ParamList;
Pnull.pts = [];
Pnull.epsi = [];

end
