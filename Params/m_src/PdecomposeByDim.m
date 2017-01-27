function Ps = PdecomposeByDim(P, dim)
% PdecomposeByDim returns an array of params sets in which P.pts(dim,:) is
%                 unique

if iscell(dim)
  dim = FindParam(P,dim);
end

X = P.pts(dim,:);
[~,IA,IC] = unique(X','rows');

for ia = 1:numel(IA)
  idx =  find(IC==ia);
  Ps(ia) = Sselect(P,idx); 
end

