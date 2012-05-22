function S = SetParam(S,ParamList,ParamValues)
% SETPARAM set the values of parameters in a parameter set
% 
% Synopsis: ParamValues = SetParam(P, ParamList, ParamValues)
%  
% Example (for Lorenz84 system):
%  
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a' );
%    val10 = 10*val; 
%    Pr10 = SetParam(Pr, 'a', val10); % values for 'a' in Pr10 are ten times those in Pr
%  
% SEE ALSO GETPARAM
%  
 
  if isstr(ParamList) 
    ind = FindParam(S,ParamList);
    S.pts(ind,:) = ParamValues;
    return
  end
  
  if isnumeric(ParamList)
    for i= 1:numel(ParamList)
      S.pts(ParamList(i),:) = ParamValues(i,:);
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
  
