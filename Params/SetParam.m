function S = SetParam(S,ParamList,ParamValues)

  if isstr(ParamList) 
    ind = FindParam(S,ParamList);
    S.pts(ind,:) = ParamValues;
    return
  end
  
  if isnumeric(ParamList)
    for i= 1:numel(ParamList)
      S.pts(ParamList(i),:) = ParamValues(i,:)
    end
    return
  end
  
  for i= 1:numel(ParamList)
    ind = FindParam(S,ParamList{i});
    S.pts(ind,:) =  ParamValues(i,:);    
  end
  
function index=  FindParam(Sys,param)
  
  if ~isfield(Sys,'ParamList')
    error('No parameter list ...');
  end

  for j = 1:numel(Sys.ParamList)
    if (strcmp(Sys.ParamList{j}, param))
      index = j;
      return;            
    end
  end
  
  error(['Parameter ' param ' not found']);
  
