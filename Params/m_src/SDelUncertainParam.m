function S = SDelUncertainParam(S, is,dims)
  
% does not remove dims 
  
  if ~exist('dims')
    dims = [];
  end
    
  
  for (i = is)
    if (isempty(find(dims==i)))
      idx = find(S.dim~=i);
      S.dim = S.dim(idx);
      S.epsi = S.epsi(idx,:);    
    end
  end
  