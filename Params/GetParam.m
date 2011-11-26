function ParamValues = GetParam(S,ParamList)
% GETPARAM get the values of parameters in parameter set
%

  if (isfield(S,'pts'))
    p= S.pts(:,1);
  else
    p = S.p;
  end
  
  if isstr(ParamList) 
    ind = FindParam(S,ParamList);    
    ParamValues = p(ind,:);
    return
  end
  
  if isnumeric(ParamList)
    ParamValues = [];  
    for i= 1:numel(ParamList)
      ParamValues(i,:) = p(ParamList(i),:);    
    end
    return
  end
    
  ParamValues = [];  
  for i= 1:numel(ParamList)
    ind = FindParam(S,ParamList{i});
    ParamValues(i,:) = p(ind,:);    
   end
    
  
